package de.charite.compbio.jannovar.hgvs.parser;

import org.antlr.v4.runtime.ANTLRInputStream;
import org.antlr.v4.runtime.BaseErrorListener;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.RecognitionException;
import org.antlr.v4.runtime.Recognizer;
import org.antlr.v4.runtime.Token;
import org.antlr.v4.runtime.tree.ParseTree;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import de.charite.compbio.jannovar.hgvs.HGVSVariant;
import de.charite.compbio.jannovar.hgvs.legacy.LegacyVariant;

/**
 * Parser for legacy mutations syntax (starting with "IVS", "EX", or "E").
 * 
 * @author Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>
 */
public class LegacyMutationParser {

	private static final Logger LOGGER = LoggerFactory.getLogger(LegacyMutationParser.class);

	private boolean debug = false;

	public LegacyMutationParser() {
	}

	public LegacyMutationParser(boolean debug) {
		this.debug = debug;
	}

	public LegacyVariant parseLegacyMutationString(String inputString) {
		LOGGER.trace("Parsing input string " + inputString);
		Antlr4HGVSParser parser = getParser(inputString);
		Antlr4HGVSParserListenerImpl listener = new Antlr4HGVSParserListenerImpl();
		parser.addParseListener(listener);
		parser.setTrace(true);
		ParseTree tree = parser.legacy_variant();
		if (debug)
			System.err.println(tree.toStringTree(parser));
		return listener.getLegacyVariant();
	}

	private Antlr4HGVSParser getParser(String inputString) {
		if (debug) {
			ANTLRInputStream inputStream = new ANTLRInputStream(inputString);
			Antlr4HGVSLexer l = new Antlr4HGVSLexer(inputStream);
			System.err.println(l.getAllTokens());
		}
		if (debug) {
			Antlr4HGVSLexer lexer = new Antlr4HGVSLexer(new ANTLRInputStream(inputString));
			// lexer.pushMode(mode);
			System.err.println("Lexer tokens");
			for (Token t : lexer.getAllTokens())
				System.err.println("\t" + t.getText() + "\t" + t);
			System.err.println("END OF LEXER TOKENS");
		}
		ANTLRInputStream inputStream = new ANTLRInputStream(inputString);
		Antlr4HGVSLexer l = new Antlr4HGVSLexer(inputStream);
		// l.pushMode(mode);
		Antlr4HGVSParser p = new Antlr4HGVSParser(new CommonTokenStream(l));
		p.setTrace(debug);
		p.addErrorListener(new BaseErrorListener() {
			@Override
			public void syntaxError(Recognizer<?, ?> recognizer, Object offendingSymbol, int line,
					int charPositionInLine, String msg, RecognitionException e) {
				throw new IllegalStateException("failed to parse at line " + line + " due to " + msg, e);
			}
		});
		return p;
	}

}