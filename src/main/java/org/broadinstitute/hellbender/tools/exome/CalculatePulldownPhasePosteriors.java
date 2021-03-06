package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.tools.exome.allelefraction.*;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountWithPhasePosteriors;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountWithPhasePosteriorsCollection;

import org.broadinstitute.hellbender.tools.exome.pulldown.Pulldown;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.ParameterReader;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Given a {@link ACNVModeledSegment} file and a file containing a {@link PosteriorSummary} for each global parameter
 * of the {@link AlleleFractionState} (both generated by fitting {@link AlleleFractionData} with an
 * {@link AlleleFractionModeller}), outputs a {@link Pulldown} that gives the probability
 * for each het to be ref minor, alt minor, or an outlier, according to the MAP model fit.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "(EXPERIMENTAL) Given case-sample output generated by the GATK ACNV workflow, outputs a pulldown file giving " +
                "ref/alt counts at heterozygous SNP sites, along with the probabilities of each site being ref minor, " +
                "alt minor, or an outlier, according to the allele-fraction model fit by ACNV.",
        oneLineSummary = "(EXPERIMENTAL) Calculate a pulldown with phase posteriors using allelic-count data and GATK ACNV output ",
        programGroup = CopyNumberProgramGroup.class
)
@BetaFeature
public class CalculatePulldownPhasePosteriors extends CommandLineProgram {
    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage or GetBayesianHetCoverage tools).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Input file for GATK ACNV segments.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File segmentsFile;

    @Argument(
            doc = "Input file for GATK ACNV allele-fraction global parameters " +
                    "(with extension " + AllelicCNV.AF_PARAMETER_FILE_SUFFIX + ").",
            fullName = ExomeStandardArgumentDefinitions.AF_PARAMETER_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.AF_PARAMETER_FILE_SHORT_NAME,
            optional = false
    )
    protected File parametersFile;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPoNFile;

    @Argument(
            doc = "Output file.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    protected File outputFile;

    @Override
    public Object doWork() {
        if (!new HDF5Library().load(null)) {  //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }

        //read counts, segments, and parameters from files
        final AllelicCountCollection counts = new AllelicCountCollection(snpCountsFile);
        final List<ACNVModeledSegment> segments = SegmentUtils.readACNVModeledSegmentFile(segmentsFile);
        final AlleleFractionState state = reconstructState(segments, parametersFile);

        //load allelic-bias panel of normals if provided
        final AllelicPanelOfNormals allelicPoN =
                allelicPoNFile != null ? AllelicPanelOfNormals.read(allelicPoNFile) : AllelicPanelOfNormals.EMPTY_PON;

        //calculate phase posteriors
        final List<SimpleInterval> unmodeledSegments = segments.stream().map(ACNVModeledSegment::getInterval).collect(Collectors.toList());
        final AllelicCountWithPhasePosteriorsCollection countsWithPhasePosteriors =
                calculatePhasePosteriors(counts, unmodeledSegments, state, allelicPoN);

        //write phase posteriors to file with same verbosity as input file
        countsWithPhasePosteriors.write(outputFile, counts.getVerbosity());

        return "SUCCESS";
    }

    private AlleleFractionState reconstructState(final List<ACNVModeledSegment> segments, final File parametersFile) {
        final Map<AlleleFractionParameter, PosteriorSummary> parameterMap = new LinkedHashMap<>();
        try (final ParameterReader<AlleleFractionParameter> reader = new ParameterReader<>(parametersFile, AlleleFractionParameter.class)) {
            reader.stream().forEach(p -> parameterMap.put(p.getKey(), p.getValue()));
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(parametersFile);
        }
        final double meanBias = parameterMap.get(AlleleFractionParameter.MEAN_BIAS).getCenter();
        final double biasVariance = parameterMap.get(AlleleFractionParameter.BIAS_VARIANCE).getCenter();
        final double outlierProbability = parameterMap.get(AlleleFractionParameter.OUTLIER_PROBABILITY).getCenter();
        final List<Double> posteriorModeMinorFractions =
                segments.stream().map(s -> s.getMinorAlleleFractionPosteriorSummary().getCenter()).collect(Collectors.toList());
        final AlleleFractionState.MinorFractions minorFractions = new AlleleFractionState.MinorFractions(posteriorModeMinorFractions);
        return new AlleleFractionState(meanBias, biasVariance, outlierProbability, minorFractions);
    }

    @VisibleForTesting
    protected static AllelicCountWithPhasePosteriorsCollection calculatePhasePosteriors(final AllelicCountCollection counts,
                                                                                        final List<SimpleInterval> segments,
                                                                                        final AlleleFractionState state,
                                                                                        final AllelicPanelOfNormals allelicPoN) {
        final TargetCollection<SimpleInterval> segmentTargetCollection = new HashedListTargetCollection<>(segments);
        final AllelicCountWithPhasePosteriorsCollection countsWithPhasePosteriors = new AllelicCountWithPhasePosteriorsCollection();
        for (final AllelicCount count : counts.getCounts()) {
            final int segmentIndex = segmentTargetCollection.index(count.getInterval());
            if (segmentIndex < 0) {
                throw new UserException.EmptyIntersection(String.format("The AllelicCount at %s is not located within one of the input segments.", count.getInterval()));
            }
            final AlleleFractionGlobalParameters parameters = state.globalParameters();
            final double minorFraction = state.segmentMinorFraction(segmentIndex);
            final double refMinorLogProb = AlleleFractionLikelihoods.hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.REF_MINOR, allelicPoN);
            final double altMinorLogProb = AlleleFractionLikelihoods.hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.ALT_MINOR, allelicPoN);
            final double outlierLogProb = AlleleFractionLikelihoods.hetLogLikelihood(parameters, minorFraction, count, AlleleFractionIndicator.OUTLIER, allelicPoN);
            final AllelicCountWithPhasePosteriors countWithPhasePosteriors = new AllelicCountWithPhasePosteriors(count, refMinorLogProb, altMinorLogProb, outlierLogProb);
            countsWithPhasePosteriors.add(countWithPhasePosteriors);
        }
        return countsWithPhasePosteriors;
    }
}
