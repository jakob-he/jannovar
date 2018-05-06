package de.charite.compbio.jannovar.cmd.hgvs_to_vcf;

import java.io.*;
import java.io.BufferedReader;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.net.ServerSocket;
import java.net.Socket;

import com.google.common.collect.Lists;

import de.charite.compbio.jannovar.JannovarException;
import de.charite.compbio.jannovar.UncheckedJannovarException;
import de.charite.compbio.jannovar.cmd.CommandLineParsingException;
import de.charite.compbio.jannovar.cmd.JannovarAnnotationCommand;
import de.charite.compbio.jannovar.hgvs.HGVSVariant;
import de.charite.compbio.jannovar.hgvs.bridge.CannotTranslateHGVSVariant;
import de.charite.compbio.jannovar.hgvs.bridge.NucleotideChangeToGenomeVariantTranslator;
import de.charite.compbio.jannovar.hgvs.nts.variant.SingleAlleleNucleotideVariant;
import de.charite.compbio.jannovar.hgvs.parser.HGVSParser;
import de.charite.compbio.jannovar.hgvs.parser.HGVSParsingException;
import de.charite.compbio.jannovar.reference.GenomeVariant;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.core.jmx.Server;

/**
 * Project transcript to chromosomal changes
 *
 * @author <a href="mailto:manuel.holtgrewe@bihealth.de">Manuel Holtgrewe</a>
 */
public class ProjectTranscriptToChromosome extends JannovarAnnotationCommand {

	/** Configuration */
	private ProjectTranscriptToChromosomeOptions options;

	/** FAI-indexed FASTA file to use */
	IndexedFastaSequenceFile fasta;

	/** Translation of variants */
	NucleotideChangeToGenomeVariantTranslator translator;

	public ProjectTranscriptToChromosome(Namespace args) throws CommandLineParsingException {
		this.options = new ProjectTranscriptToChromosomeOptions();
		this.options.setFromArgs(args);
	}

	@Override
	public void run() throws JannovarException {
		System.err.println("Options");
		System.err.println(options.toString());
		System.err.println("Loading database file...");
		deserializeTranscriptDefinitionFile(options.getDatabaseFilePath());
		System.err.println("Loading FASTA index...");
		loadFASTAIndex();
		if (options.getServerMode()) {
			serverProcessing();
		} else {
			simpleFileProcessing();
		}
	}

	private void serverProcessing() {
		System.err.println("Running HGVS to VCF translator in server mode.");
		System.err.println("Running on port: " + Integer.toString(options.getServerPort()));
		try (ServerSocket server = new ServerSocket(options.getServerPort())) {
			while (true) {
				acceptConnection(server);
			}
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}

	private void acceptConnection(ServerSocket server) {
		try (
			Socket client = server.accept();
			BufferedReader in = new BufferedReader(new InputStreamReader(client.getInputStream()));
			PrintWriter pw = new PrintWriter(client.getOutputStream());
		) {
			/* get message size in first line */
			String b;
			int read, messageSize = 0, recvSize = 0;
			if ((b = in.readLine()) != null) {
				messageSize = Integer.parseInt(b);
			}
			char[] chunk = new char[2048];
			StringBuilder chunks = new StringBuilder();
			while (recvSize < messageSize) {
				if ((read = in.read(chunk, recvSize, Integer.min(messageSize - recvSize, 2048))) < 1) {
					break;
				}
				recvSize += read;
				chunks.append(chunk);
			}

			if (messageSize > 0) {
				ByteArrayOutputStream output = new ByteArrayOutputStream();
				BufferedReader bsr = new BufferedReader(new StringReader(chunks.toString()));
				try (VariantContextWriter writer = openStream(output)) {
					processVariants(writer, bsr);
					pw.println(output.size());
					pw.println(0);
					pw.print(output.toString());
				} catch (Exception e) {
					String msg = e.getMessage();
					pw.println(msg.length() + 1);
					pw.println(-1);
					pw.println(e.getMessage());
				}
			} else {
				System.err.println("No message received.");
			}

			System.err.println("Connection closed.");
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}

	private void simpleFileProcessing(){
		System.err.println("Opening output VCF file...");
		try (VariantContextWriter writer = openOutputFile()) {
			processFile(writer);
		}
	}

	private VariantContextWriter openStream(OutputStream output) {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
		builder.setReferenceDictionary((fasta.getSequenceDictionary()));
		builder.setOutputStream(output);
		builder.unsetOption(Options.INDEX_ON_THE_FLY);
		VariantContextWriter writer = builder.build();
		VCFHeader header = createVCFHeader();
		writer.writeHeader(header);
		return writer;
	}

	private VariantContextWriter openOutputFile() {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
			.setReferenceDictionary(fasta.getSequenceDictionary()).setOutputFile(options.getPathOutputVCF());
		if (options.getPathOutputVCF().endsWith(".gz") || options.getPathOutputVCF().endsWith(".bcf"))
			builder.setOption(Options.INDEX_ON_THE_FLY);
		else
			builder.unsetOption(Options.INDEX_ON_THE_FLY);
		VariantContextWriter writer = builder.build();

		VCFHeader header = createVCFHeader();

		writer.writeHeader(header);

		return writer;
	}

	private VCFHeader createVCFHeader() {
		VCFHeader header = new VCFHeader();
		int i = 0;
		for (SAMSequenceRecord record : fasta.getSequenceDictionary().getSequences()) {
			Map<String, String> mapping = new TreeMap<String, String>();
			mapping.put("ID", record.getSequenceName());
			mapping.put("length", Integer.toString(record.getSequenceLength()));
			header.addMetaDataLine(new VCFContigHeaderLine(mapping, i++));
		}

		header.addMetaDataLine(new VCFSimpleHeaderLine("ALT", "ERROR", "Error in conversion"));
		header.addMetaDataLine(new VCFFilterHeaderLine("PARSE_ERROR",
			"Problem in parsing original HGVS variant string, written out as variant at 1:g.1N>N"));
		header.addMetaDataLine(new VCFInfoHeaderLine("ERROR_MESSAGE", 1, VCFHeaderLineType.String, "Error message"));
		header.addMetaDataLine(new VCFInfoHeaderLine("ORIG_VAR", 1, VCFHeaderLineType.String,
			"Original HGVS variant string from input file to hgvs-to-vcf"));
		return header;
	}

	private void loadFASTAIndex() {
		try {
			this.fasta = new IndexedFastaSequenceFile(new File(options.getPathReferenceFASTA()));
		} catch (FileNotFoundException e) {
			throw new UncheckedJannovarException("Could not load FASTA index", e);
		}
		if (this.fasta.getSequenceDictionary() == null) {
			throw new UncheckedJannovarException(
					"FASTA sequence dictionary empty, you have a REFERENCE.dict file (create with Picard "
							+ "or samtools dict, version >=1.3)");
		}

		this.translator = new NucleotideChangeToGenomeVariantTranslator(jannovarData, fasta);
	}

	/**
	 * @return Variant context indicating error
	 */
	private VariantContext buildErrorVariantContext(String origString, String message) {
		Allele alleleRef = Allele.create("N", true);
		Allele alleleAlt = Allele.create("<ERROR>", false);
		return new VariantContextBuilder().loc("1", 1, 1).alleles(Lists.newArrayList(alleleRef, alleleAlt))
				.filter("PARSE_ERROR").attribute("ORIG_VAR", urlEncode(origString))
				.attribute("ERROR_MESSAGE", urlEncode(message)).make();
	}

	private String urlEncode(String s) {
		try {
			return URLEncoder.encode(s, "utf-8").replaceAll("=", "%3D");
		} catch (UnsupportedEncodingException e) {
			return s;
		}
	}

	private void processVariants(VariantContextWriter writer, BufferedReader br) {
		String line;
		try {
			while ((line = br.readLine()) != null) {
				// Read line
				String word = line.trim();

				// stop if empty line encountered
				if (word.isEmpty()) {
					break;
				}
				// Parse variant
				HGVSVariant rawVar;
				try {
					HGVSParser parser = new HGVSParser();
					rawVar = parser.parseHGVSString(word);
					if (!(rawVar instanceof SingleAlleleNucleotideVariant)) {
						writer.add(buildErrorVariantContext(word, "More than one allele in nucleotide variant"));
						continue;
					}
				} catch (HGVSParsingException e) {
					writer.add(buildErrorVariantContext(word, e.getMessage()));
					continue;
				}

				// Convert from transcript to genome variant
				GenomeVariant genomeVar = translate((SingleAlleleNucleotideVariant) rawVar);
				if (genomeVar == null) {
					writer.add(buildErrorVariantContext(word, "Could not translate HGVS to genomic variant"));
					continue;
				}

				// Write variant to VCF file
				writeVariant(writer, genomeVar);
				System.err.println(word + " => " + rawVar + " => " + genomeVar);
			}
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}

	private void processFile(VariantContextWriter writer) {
		try (BufferedReader br = new BufferedReader(new FileReader(new File(options.getPathInputText())))) {
			processVariants(writer, br);
		} catch (FileNotFoundException e) {
			throw new UncheckedJannovarException("Problem opening file", e);
		} catch (IOException e) {
			throw new UncheckedJannovarException("Problem reading from file", e);
		}
	}

	/** Map contig name (from genome variant) to contig name in FASTA */
	private String mapContigToFasta(String contigName) {
		// Map genome variant's contig to unique ID
		Integer contigID = jannovarData.getRefDict().getContigNameToID().get(contigName);
		if (contigID == null)
			throw new UncheckedJannovarException("Unknown contig name " + contigName);
		// Try to find matching contig in fasta
		String nameInFasta = null;
		for (SAMSequenceRecord record : fasta.getSequenceDictionary().getSequences()) {
			if (jannovarData.getRefDict().getContigNameToID().containsKey(record.getSequenceName())) {
				nameInFasta = record.getSequenceName();
				break;
			}
		}
		if (nameInFasta == null)
			throw new UncheckedJannovarException("Could not find corresponding contig in FASTA for " + contigName);

		return nameInFasta;
	}

	private void writeVariant(VariantContextWriter writer, GenomeVariant genomeVar) {
		String nameInFasta = mapContigToFasta(genomeVar.getChrName());
		List<Allele> alleles = new ArrayList<Allele>();
		int shift = 0;
		if (genomeVar.getRef().isEmpty() || genomeVar.getAlt().isEmpty()) {
			shift = -1;
			String left = fasta.getSubsequenceAt(nameInFasta, genomeVar.getPos(), genomeVar.getPos())
					.getBaseString();
			alleles.add(Allele.create(left + genomeVar.getRef(), true));
			alleles.add(Allele.create(left + genomeVar.getAlt(), false));
		} else {
			alleles.add(Allele.create(genomeVar.getRef(), true));
			alleles.add(Allele.create(genomeVar.getAlt(), false));
		}

		VariantContextBuilder builder = new VariantContextBuilder();
		builder.chr(genomeVar.getChrName()).start(genomeVar.getPos() + shift + 1)
				.computeEndFromAlleles(alleles, genomeVar.getPos() + shift + 1).alleles(alleles);

		writer.add(builder.make());
	}

	private GenomeVariant translate(SingleAlleleNucleotideVariant rawVar) {
		try {
			return translator.translateNucleotideVariantToGenomeVariant((SingleAlleleNucleotideVariant) rawVar, true);
		} catch (CannotTranslateHGVSVariant e) {
			System.err.println("Could not translate variant " + rawVar + ": " + e.toString());
			return null;
		}
	}

}
