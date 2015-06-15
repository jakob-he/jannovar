package de.charite.compbio.jannovar.hgvs.parser.nts.change;

import org.junit.Assert;
import org.junit.Test;

import de.charite.compbio.jannovar.hgvs.parser.HGVSLexer;
import de.charite.compbio.jannovar.hgvs.parser.HGVSParser;
import de.charite.compbio.jannovar.hgvs.parser.HGVSParserTestBase;
import de.charite.compbio.jannovar.hgvs.parser.HGVSParser.Nt_change_substitutionContext;

/**
 * Parser for HGVS substitution nucleotide changes.
 * 
 * @author Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>
 */
public class HGVSParserNucleotideMiscChangeTest extends HGVSParserTestBase {

	@Test
	public void testNucleotideSubstitution() {
		HGVSParser parser = buildParserForString("123C>A", HGVSLexer.NUCLEOTIDE_CHANGE, false);
		Nt_change_substitutionContext nt_change_substitution = parser.nt_change_substitution();
		Assert.assertEquals("(nt_change_substitution (nt_point_location (nt_base_location (nt_number 123))) C > A)",
				nt_change_substitution.toStringTree(parser));
	}

}