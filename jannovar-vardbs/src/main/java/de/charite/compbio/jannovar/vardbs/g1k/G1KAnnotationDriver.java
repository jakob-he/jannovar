package de.charite.compbio.jannovar.vardbs.g1k;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import de.charite.compbio.jannovar.vardbs.base.*;
import de.charite.compbio.jannovar.vardbs.base.vcf.AbstractVCFDBAnnotationDriver;
import de.charite.compbio.jannovar.vardbs.base.vcf.GenotypeMatch;
import de.charite.compbio.jannovar.vardbs.dbsnp.DBSNPRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;

// TODO: handle MNVs appropriately

/**
 * Annotation driver class for annotations using 1KG
 *
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 */
public class G1KAnnotationDriver extends AbstractVCFDBAnnotationDriver<G1KRecord> {

	public G1KAnnotationDriver(String vcfPath, String fastaPath, DBAnnotationOptions options)
			throws JannovarVarDBException {
		super(new VCFReaderVariantProvider(vcfPath), fastaPath, options, new G1KVariantContextToRecordConverter());
	}

	@Override
	public VCFHeaderExtender constructVCFHeaderExtender() {
		return new G1KVCFHeaderExtender(options);
	}

	@Override
	protected HashMap<Integer, AnnotatingRecord<G1KRecord>> pickAnnotatingDBRecords(
			HashMap<Integer, ArrayList<GenotypeMatch>> annotatingRecords,
			HashMap<GenotypeMatch, AnnotatingRecord<G1KRecord>> matchToRecord, boolean isMatch) {
		// Pick best annotation for each alternative allele
		HashMap<Integer, AnnotatingRecord<G1KRecord>> annotatingExacRecord = new HashMap<>();
		for (Entry<Integer, ArrayList<GenotypeMatch>> entry : annotatingRecords.entrySet()) {
			final int alleleNo = entry.getKey();
			for (GenotypeMatch m : entry.getValue()) {
				if (!annotatingExacRecord.containsKey(alleleNo)) {
					annotatingExacRecord.put(alleleNo, matchToRecord.get(m));
				} else {
					final G1KRecord current = annotatingExacRecord.get(alleleNo).getRecord();
					final G1KRecord update = matchToRecord.get(m).getRecord();
					if (update.getAlleleFrequencies(G1KPopulation.ALL).size() < alleleNo)
						continue;
					else if (isMatch && current.getAlleleFrequencies(G1KPopulation.ALL).size() < alleleNo
							|| isMatch && current.highestAlleleFreq(alleleNo - 1) < update.highestAlleleFreq(alleleNo - 1))
						annotatingExacRecord.put(alleleNo, matchToRecord.get(m));
				}
			}
		}
		return annotatingExacRecord;
	}

	@Override
	protected VariantContext annotateWithDBRecords(VariantContext vc,
			HashMap<Integer, AnnotatingRecord<G1KRecord>> matchRecords,
			HashMap<Integer, AnnotatingRecord<G1KRecord>> overlapRecords) {
		if (matchRecords.isEmpty())
			return vc;

		VariantContextBuilder builder = new VariantContextBuilder(vc);

		// Annotate with records with matching allele
		annotateAlleleCounts(vc, "", matchRecords, builder);
		annotateChromosomeCounts(vc, "", matchRecords, builder);
		annotateFrequencies(vc, "", matchRecords, builder);
		annotateBestAF(vc, "", matchRecords, builder);

		// Annotate with records with overlapping positions
		if (options.isReportOverlapping() && !options.isReportOverlappingAsMatching()) {
			annotateAlleleCounts(vc, "OVL_", overlapRecords, builder);
			annotateChromosomeCounts(vc, "OVL_", overlapRecords, builder);
			annotateFrequencies(vc, "OVL_", overlapRecords, builder);
			annotateBestAF(vc, "OVL_", overlapRecords, builder);
		}

		return builder.make();
	}

	private void annotateBestAF(VariantContext vc, String infix, HashMap<Integer, AnnotatingRecord<G1KRecord>> records,
			VariantContextBuilder builder) {
		ArrayList<Double> afs = new ArrayList<>();
		ArrayList<Integer> acs = new ArrayList<>();
		for (int i = 1; i < vc.getNAlleles(); ++i) {
			if (records.get(i) == null) {
				afs.add(0.0);
				acs.add(0);
			} else {
				final G1KRecord record = records.get(i).getRecord();
				final int alleleNo = records.get(i).getAlleleNo();
				final G1KPopulation pop = record.popWithHighestAlleleFreq(alleleNo - 1);
				afs.add(record.getAlleleFrequencies(pop).get(alleleNo - 1));
				acs.add(record.getAlleleCounts().get(alleleNo - 1));
			}
		}

		builder.attribute(options.getVCFIdentifierPrefix() + infix + "BEST_AF", afs);
	}

	private void annotateChromosomeCounts(VariantContext vc, String infix,
			HashMap<Integer, AnnotatingRecord<G1KRecord>> records, VariantContextBuilder builder) {
		if (records.isEmpty())
			return;
		G1KRecord first = records.values().iterator().next().getRecord();
		builder.attribute(options.getVCFIdentifierPrefix() + infix + "AN", first.getChromCounts());
	}

	private void annotateAlleleCounts(VariantContext vc, String infix,
			HashMap<Integer, AnnotatingRecord<G1KRecord>> records, VariantContextBuilder builder) {
		final String attrID = options.getVCFIdentifierPrefix() + infix + "AC";
		ArrayList<Integer> acList = new ArrayList<>();
		for (int i = 1; i < vc.getNAlleles(); ++i) {
			if (records.get(i) == null) {
				acList.add(0);
				continue;
			}
			final G1KRecord record = records.get(i).getRecord();
			final int alleleNo = records.get(i).getAlleleNo();
			if (record.getAlleleCounts().isEmpty()) {
				acList.add(0);
			} else {
				acList.add(record.getAlleleCounts().get(alleleNo - 1));
			}

			if (!acList.isEmpty())
				builder.attribute(attrID, acList);

			builder.attribute(attrID, acList);
		}
	}

	private void annotateFrequencies(VariantContext vc, String infix,
			HashMap<Integer, AnnotatingRecord<G1KRecord>> records, VariantContextBuilder builder) {
		for (G1KPopulation pop : G1KPopulation.values()) {
			final String attrID = options.getVCFIdentifierPrefix() + infix + "AF_" + pop;
			ArrayList<Double> afList = new ArrayList<>();
			for (int i = 1; i < vc.getNAlleles(); ++i) {
				if (records.get(i) == null) {
					afList.add(0.0);
					continue;
				}
				final G1KRecord record = records.get(i).getRecord();
				final int alleleNo = records.get(i).getAlleleNo();
				if (record.getAlleleCounts().isEmpty()) {
					afList.add(0.0);
				} else {
					afList.add(record.getAlleleFrequencies(pop).get(alleleNo - 1));
				}
			}

			if (!afList.isEmpty())
				builder.attribute(attrID, afList);

			builder.attribute(attrID, afList);
		}
	}

}
