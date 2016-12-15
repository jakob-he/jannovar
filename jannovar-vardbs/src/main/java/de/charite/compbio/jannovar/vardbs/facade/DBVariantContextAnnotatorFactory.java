package de.charite.compbio.jannovar.vardbs.facade;

import de.charite.compbio.jannovar.vardbs.base.DBAnnotationOptions;
import de.charite.compbio.jannovar.vardbs.base.JannovarVarDBException;
import de.charite.compbio.jannovar.vardbs.dbsnp.DBSNPAnnotationDriver;
import de.charite.compbio.jannovar.vardbs.exac.ExacAnnotationDriver;
import de.charite.compbio.jannovar.vardbs.g1k.G1KAnnotationDriver;
import de.charite.compbio.jannovar.vardbs.remm.ReMMAnnotationDriver;
import de.charite.compbio.jannovar.vardbs.tabix.TabixAnnotationDriver;
import de.charite.compbio.jannovar.vardbs.uk10k.UK10KAnnotationDriver;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Factory for generating {@link DBVariantContextAnnotator} objects
 * 
 * @author <a href="mailto:manuel.holtgrewe@bihealth.de">Manuel Holtgrewe</a>
 */
public class DBVariantContextAnnotatorFactory {

	/**
	 * Construct dbSNP {@link VariantContext} annotator factory.
	 * 
	 * @param vcfDBPath
	 *            Path to gzip-compressed, normalized and tbi-indexed dbSNP file to use for the annotation
	 * @param fastaRefPath
	 *            Path to reference FASTA file
	 * @param options
	 *            Configuration for the variant context annotaiton
	 * @return Preconfigured {@link DBVariantContextAnnotator} object
	 * @throws JannovarVarDBException
	 *             on problems loading the resources
	 */
	public DBVariantContextAnnotator constructDBSNP(String vcfDBPath, String fastaRefPath, DBAnnotationOptions options)
			throws JannovarVarDBException {
		return new DBVariantContextAnnotator(new DBSNPAnnotationDriver(vcfDBPath, fastaRefPath, options), options);
	}

	/**
	 * Construct ExAC {@link VariantContext} annotator factory.
	 * 
	 * @param vcfDBPath
	 *            Path to gzip-compressed, normalized and tbi-indexed ExAC file to use for the annotation
	 * @param fastaRefPath
	 *            Path to reference FASTA file
	 * @param options
	 *            Configuration for the variant context annotaiton
	 * @return Preconfigured {@link DBVariantContextAnnotator} object
	 * @throws JannovarVarDBException
	 *             on problems loading the resources
	 */
	public DBVariantContextAnnotator constructExac(String vcfDBPath, String fastaRefPath, DBAnnotationOptions options)
			throws JannovarVarDBException {
		return new DBVariantContextAnnotator(new ExacAnnotationDriver(vcfDBPath, fastaRefPath, options), options);
	}

	/**
	 * Construct UK10K {@link VariantContext} annotator factory.
	 * 
	 * @param vcfDBPath
	 *            Path to gzip-compressed, normalized and tbi-indexed UK10K file to use for the annotation
	 * @param fastaRefPath
	 *            Path to reference FASTA file
	 * @param options
	 *            Configuration for the variant context annotaiton
	 * @return Preconfigured {@link DBVariantContextAnnotator} object
	 * @throws JannovarVarDBException
	 *             on problems loading the resources
	 */
	public DBVariantContextAnnotator constructUK10K(String vcfDBPath, String fastaRefPath, DBAnnotationOptions options)
			throws JannovarVarDBException {
		return new DBVariantContextAnnotator(new UK10KAnnotationDriver(vcfDBPath, fastaRefPath, options), options);
	}

	/**
	 * Construct 1KG {@link VariantContext} annotator factory.
	 * 
	 * @param vcfDBPath
	 *            Path to gzip-compressed, normalized and tbi-indexed G1K file to use for the annotation
	 * @param fastaRefPath
	 *            Path to reference FASTA file
	 * @param options
	 *            Configuration for the variant context annotaiton
	 * @return Preconfigured {@link DBVariantContextAnnotator} object
	 * @throws JannovarVarDBException
	 *             on problems loading the resources
	 */
	public DBVariantContextAnnotator construct1KG(String vcfDBPath, String fastaRefPath, DBAnnotationOptions options)
			throws JannovarVarDBException {
		return new DBVariantContextAnnotator(new G1KAnnotationDriver(vcfDBPath, fastaRefPath, options), options);
	}

	/**
	 * Construct Tabix {@link VariantContext} annotator factory.
	 * 
	 * @param tabixDBPath
	 *            Path to gzip-compressed, tbi-indexed tabix including ref and alt
	 * @param fastaRefPath
	 *            Path to reference FASTA file
	 * @param options
	 *            Configuration for the variant context annotaiton
	 * @return Preconfigured {@link DBVariantContextAnnotator} object
	 * @throws JannovarVarDBException
	 *             on problems loading the resources
	 */
	public DBVariantContextAnnotator constructTabix(String tabixDBPath, String fastaRefPath,
			DBAnnotationOptions options) throws JannovarVarDBException {
		return new DBVariantContextAnnotator(new TabixAnnotationDriver(tabixDBPath, fastaRefPath, options), options);
	}

	/**
	 * Construct ReMM {@link VariantContext} annotator factory.
	 * 
	 * @param remmDBPath
	 *            Path to gzip-compressed, tbi-indexed tabix including ref and alt
	 * @param fastaRefPath
	 *            Path to reference FASTA file
	 * @param options
	 *            Configuration for the variant context annotaiton
	 * @return Preconfigured {@link DBVariantContextAnnotator} object
	 * @throws JannovarVarDBException
	 *             on problems loading the resources
	 */
	public DBVariantContextAnnotator constructReMM(String remmDBPath, DBAnnotationOptions options)
			throws JannovarVarDBException {
		return new DBVariantContextAnnotator(new ReMMAnnotationDriver(remmDBPath, options), options);
	}

}
