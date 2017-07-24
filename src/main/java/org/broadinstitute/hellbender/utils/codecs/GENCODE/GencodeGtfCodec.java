package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jonn on 7/21/17.
 */
final public class GencodeGtfCodec extends AsciiFeatureCodec<GencodeGtfFeature> {

    private static final String COMMENT_START = "##";

    private static final int NUM_COLUMNS = 9;

    protected long currentLineNum = 1;
    protected List<String> header = new ArrayList<>();
    private final ArrayList< String[] > intermediateRecordStore = new ArrayList<>();

    public GencodeGtfCodec() {
        super(GencodeGtfFeature.class);
    }


    @Override
    public GencodeGtfFeature decode(String s) {

        GencodeGtfFeature decodedFeature = null;

        String[] splitLine = s.split("\t");

        // Ensure the file is at least trivially well-formed:
        if (splitLine.length != NUM_COLUMNS) {
            throw new UserException.MalformedFile("Found an invalid number of columns in the given GENCODE file on line "
                    + currentLineNum + " - Given: " + splitLine.length + " Expected: " + NUM_COLUMNS);
        }

        // We need to key off the feature type to collapse our accumulated records:
        final String featureType = splitLine[2];

        // Once we see another gene or transcript, we take all accumulated records and combine them into a
        // GencodeGtfFeature.
        if (featureType.equals("gene") || featureType.equals("transcript")) {
            if ( intermediateRecordStore.size() != 0 ) {

                // OK, we go through the record and

            }
        }

        ++currentLineNum;

        return decodedFeature;
    }

    @Override
    public Object readActualHeader(LineIterator reader) {

        boolean isFirst = true;

        while ( reader.hasNext() ) {
            String line = reader.peek();

            // The file will start with commented out lines.
            // Grab them until there are no more commented out lines.
            if ( line.substring(0,2).equals(COMMENT_START) ) {
                header.add(line);
                reader.next();
                isFirst = false;
            }
            else if ( isFirst ) {
                throw new UserException.MalformedFile("GENCODE file does not have a header!");
            }
            else {
                break;
            }
        }

        return header;
    }

    @Override
    public boolean canDecode(String path) {
        return path.toLowerCase().substring(0,7).equals("gencode") && path.toLowerCase().endsWith(".gtf");
    }
}
