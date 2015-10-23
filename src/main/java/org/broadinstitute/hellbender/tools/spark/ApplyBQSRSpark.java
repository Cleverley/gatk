package org.broadinstitute.hellbender.tools.spark;

import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;

import java.io.IOException;

@CommandLineProgramProperties(summary="Apply Base Quality Recalibration to a bam file using spark",
        oneLineSummary="apply BQSR on spark",
        programGroup = SparkProgramGroup.class)
public final class ApplyBQSRSpark extends GATKSparkTool {
    private static final long serialVersionUID = 0l;

    @Override
    public boolean requiresReads() { return true; }

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    /**
     * Enables recalibration of base qualities.
     * The covariates tables are produced by the BaseRecalibrator tool.
     * Please be aware that you should only run recalibration with the covariates file created on the same input bam(s).
     */
    @Argument(fullName="bqsr_recal_file", shortName="bqsr", doc="Input covariates table file for base quality score recalibration")
    private String bqsrRecalFile;

    @ArgumentCollection
    private ApplyBQSRArgumentCollection applyBQSRArgs = new ApplyBQSRArgumentCollection();

    @Argument(doc = "If specified, shard the output bam", shortName = "shardedOutput", fullName = "shardedOutput", optional = true)
    private boolean shardedOutput = false;

    @Override
    protected void runTool(JavaSparkContext ctx) {
        JavaRDD<GATKRead> initialReads = getReads();
        Broadcast<RecalibrationReport> recalibrationReportBroadCast = ctx.broadcast(new RecalibrationReport(BucketUtils.openFile(bqsrRecalFile, getAuthHolder())));
        final JavaRDD<GATKRead> recalibratedReads = ApplyBQSRSparkFn.apply(initialReads, recalibrationReportBroadCast, getHeaderForReads(), applyBQSRArgs);

        try {
            ReadsSparkSink.writeReads(ctx, output, recalibratedReads, getHeaderForReads(), shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }
    }
}
