package org.broadinstitute.hellbender.utils.codecs.GENCODE;

import htsjdk.tribble.Feature;

/**
 * A GTF Feature represents one row of a GTF File.
 * The specification of a GTF file is defined here:
 * http://mblab.wustl.edu/GTF22.html
 *
 * Created by jonn on 7/21/17.
 */
public class GencodeGtfTranscriptFeature implements Feature {

    @Override
    public String getContig() {
        return null;
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }

    // ================================================================================================

    public enum AnnotationSource {
        ENSEMBL,
        HAVANA
    }

    public enum FeatureType {
        GENE,
        TRANSCRIPT,
        EXON,
        CDS,
        UTR,
        START_CODON,
        STOP_CODON,
        SELENOCYSTEINE
    }

    public enum GenomicStrand {
        FORWARD,
        BACKWARD
    }

    public enum GenomicPhase {
        ZERO,
        ONE,
        TWO,
        DOT
    }

    public enum GeneTranscriptType {
        // Immunoglobulin (Ig) variable chain and T-cell receptor (TcR) genes imported or annotated according to the IMGT (http://www.imgt.org/)
        IG_C_gene,
        IG_D_gene,
        IG_J_gene,
        IG_LV_gene,
        IG_V_gene,
        TR_C_gene,
        TR_J_gene,
        TR_V_gene,
        TR_D_gene,

        // Inactivated immunoglobulin gene.
        IG_pseudogene,
        IG_C_pseudogene,
        IG_J_pseudogene,
        IG_V_pseudogene,
        TR_V_pseudogene,
        TR_J_pseudogene,

        // Non-coding RNA predicted using sequences from Rfam (http://rfam.xfam.org/) and miRBase (http://www.mirbase.org/)
        Mt_rRNA,
        Mt_tRNA,
        miRNA,
        misc_RNA,
        rRNA,
        scRNA,
        snRNA,
        snoRNA,
        ribozyme,
        sRNA,
        scaRNA,

        // Non-coding RNA predicted to be pseudogene by the Ensembl pipeline
        Mt_tRNA_pseudogene,
        tRNA_pseudogene,
        snoRNA_pseudogene,
        snRNA_pseudogene,
        scRNA_pseudogene,
        rRNA_pseudogene,
        misc_RNA_pseudogene,
        miRNA_pseudogene,

        // To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies.
        TEC,

        // If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD.
        nonsense_mediated_decay,

        // Transcript that has polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation.
        non_stop_decay,

        // Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants.
        retained_intron,

        // Contains an open reading frame (ORF).
        protein_coding,

        // Doesn't contain an ORF.
        processed_transcript,

        // Transcript which is known from the literature to not be protein coding.
        non_coding,

        // Transcript believed to be protein coding, but with more than one possible open reading frame.
        ambiguous_orf,

        // Long non-coding transcript in introns of a coding gene that does not overlap any exons.
        sense_intronic,

        // Long non-coding transcript that contains a coding gene in its intron on the same strand.
        sense_overlapping,

        // Has transcripts that overlap the genomic span (i.e. exon or introns) of a protein-coding locus on the opposite strand.
        antisense,

        known_ncrna,

        // Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following.
        pseudogene,

        // Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome.
        processed_pseudogene,

        // Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated.
        polymorphic_pseudogene,

        // Pseudogene owing to a reverse transcribed and re-inserted sequence.
        retrotransposed,

        // Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression.
        transcribed_processed_pseudogene,
        transcribed_unprocessed_pseudogene,
        transcribed_unitary_pseudogene,

        // Pseudogene that has mass spec data suggesting that it is also translated.
        translated_processed_pseudogene,
        translated_unprocessed_pseudogene,

        // A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species.
        unitary_pseudogene,

        // Pseudogene that can contain introns since produced by gene duplication.
        unprocessed_pseudogene,

        // Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
        artifact,

        // Long, intervening noncoding (linc) RNA that can be found in evolutionarily conserved, intergenic regions.
        lincRNA,

        // Unspliced lncRNA that is several kb in size.
        macro_lncRNA,

        // Transcript where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR.
        three_prime_overlapping_ncRNA,

        // Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain.
        disrupted_domain,

        // Short non coding RNA gene that forms part of the vault ribonucleoprotein complex.
        vaultRNA,

        // A non-coding locus that originates from within the promoter region of a protein-coding gene, with transcription proceeding in the opposite direction on the other strand.
        bidirectional_promotor_IncRNA

    }

    public enum GeneTranscriptStatus {
        KNOWN,
        NOVEL,
        PUTATIVE
    }

    public enum Tag {
        // 3' end extended based on RNA-seq data.
        three_nested_supported_extension,

        // 3' end extended based on RNA-seq data.
        three_standard_supported_extension,

        // annotated based on RNA-seq data.
        fourfivefour_RNA_Seq_supported,

        // 5' end extended based on RNA-seq data.
        five_nested_supported_extension,

        // 5' end extended based on RNA-seq data.
        five_standard_supported_extension,

        // shares an identical CDS but has alternative 5' UTR with respect to a reference variant.
        alternative_3_UTR,

        // shares an identical CDS but has alternative 3' UTR with respect to a reference variant.
        alternative_5_UTR,

        // (This flag corresponds to the older flag "appris_principal") Where the transcript expected to code for the main
        appris_principal_1,

        // (This flag corresponds to the older flag "appris_candidate_ccds") Where the APPRIS core modules are unable to choose a
        appris_principal_2,

        // Where the APPRIS core modules are unable to choose a clear principal variant and there more than one of the variants
        appris_principal_3,

        // (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_ccds") Where the APPRIS core modules are unable
        appris_principal_4,

        // (This flag corresponds to the Ensembl 78 flag "appris_candidate_longest_seq") Where the APPRIS core modules are unable
        appris_principal_5,

        // Candidate transcript(s) models that are conserved in at least three tested non-primate species.
        appris_alternative_1,

        // Candidate transcript(s) models that appear to be conserved in fewer than three tested non-primate species.
        appris_alternative_2,

        // ranscript expected to code for the main functional isoform based on a range of protein features (APPRIS pipeline).
        appris_principal,

        // where there is no single 'appris_principal' variant the main functional isoform will be translated from one of the
        appris_candidate,

        // he "appris_candidate" transcript that has an unique CCDS.
        appris_candidate_ccds,

        // where there is no 'appris_principal' variant, the candidate with highest APPRIS score is selected as the primary
        appris_candidate_highest_score,

        // where there is no 'appris_principal' variant, the longest of the 'appris_candidate' variants is selected as the primary
        appris_candidate_longest,

        // he "appris_candidate" transcripts where there are several CCDS, in this case APPRIS labels the longest CCDS.
        appris_candidate_longest_ccds,

        // where there is no "appris_candidate_ccds" or "appris_candidate_longest_ccds" variant, the longest protein of the
        appris_candidate_longest_seq,

        // identifies a subset of representative transcripts for each gene; prioritises full-length protein coding transcripts
        basic,

        // ranscript contains two confidently annotated CDSs. Support may come from eg proteomic data, cross-species conservation
        bicistronic,

        // ranscript 5' end overlaps ENCODE or Fantom CAGE cluster.
        CAGE_supported_TSS,

        // member of the consensus CDS gene set, confirming coding regions between ENSEMBL, UCSC, NCBI and HAVANA.
        CCDS,

        // he coding region end could not be confirmed.
        cds_end_NF,

        // he coding region start could not be confirmed.
        cds_start_NF,

        // ranscript QC checked using dotplot to identify features eg splice junctions, end of homology.
        dotter_confirmed,

        // an upstream ATG is used where a downstream ATG seems more evolutionary conserved.
        downstream_ATG,

        // ranscript was tested and confirmed experimentally.
        exp_conf,

        // locus consists of non-overlapping transcript fragments either because of genome assembly issues (i.e., gaps or
        fragmented_locus,

        // ranscript model contains all possible in-frame exons supported by homology, experimental evidence or conservation, but
        inferred_exon_combination,

        // ranscript model is not supported by a single piece of transcript evidence. May be supported by multiple fragments of
        inferred_transcript_model,

        // ranscript supported by transcript evidence that, while ampping best-in-genome, shows regions of poor sequence quality.
        low_sequence_quality,

        // he mRNA end could not be confirmed.
        mRNA_end_NF,

        // he mRNA start could not be confirmed.
        mRNA_start_NF,

        // in-frame type of variation where, at the acceptor site, some variants splice after the first AG and others after the
        NAGNAG_splice_site,

        // he locus is a host for small non-coding RNAs.
        ncRNA_host,

        // annotated based on RNA-seq data.
        nested_454_RNA_Seq_supported,

        // he transcript looks like it is subject to NMD but publications, experiments or conservation support the translation of
        NMD_exception,

        // codon if the transcript were longer but cannot currently be annotated as NMD as does not fulfil all criteria - most
        NMD_likely_if_extended,

        // he CDS has a non-ATG start and its validity is supported by publication or conservation.
        non_ATG_start,

        // he transcript has a non-canonical splice site conserved in other species.
        non_canonical_conserved,

        // he transcript has a non-canonical splice site explained by a genomic sequencing error.
        non_canonical_genome_sequence_error,

        // he transcript has a non-canonical splice site explained by other reasons.
        non_canonical_other,

        // he transcript has a non-canonical splice site explained by a SNP.
        non_canonical_polymorphism,

        // he transcript has a non-canonical splice site that needs experimental confirmation.
        non_canonical_TEC,

        // he transcript has a non-canonical splice site explained by a U12 intron (i.e. AT-AC splice site).
        non_canonical_U12,

        // a splice variant for which supporting evidence has not been submitted to databases, i.e. the model is based on
        non_submitted_evidence,

        // a transcript is supported by evidence from same species paralogous loci.
        not_best_in_genome_evidence,

        // evidence from other species was used to build model.
        not_organism_supported,

        // protein-coding locus with no paralogues or orthologs.
        orphan,

        // exon(s) of the locus overlap exon(s) of a readthrough transcript or a transcript belonging to another locus.
        overlapping_locus,

        // a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high
        overlapping_uORF,

        // annotation in the pseudo-autosomal region, which is duplicated between chromosomes X and Y.
        PAR,

        // member of the pseudogene set predicted by YALE, UCSC and HAVANA.
        pseudo_consens,

        // a transcript that overlaps two or more independent loci but is considered to belong to a third, separate locus.
        readthrough_transcript,

        // locus overlaps a sequence error or an assembly error in the reference genome that affects its annotation (e.g., 1 or
        reference_genome_error,

        // internal intron of CDS portion of transcript is retained.
        retained_intron_CDS,

        // final intron of CDS portion of transcript is retained.
        retained_intron_final,

        // first intron of CDS portion of transcript is retained.
        retained_intron_first,

        // protein-coding locus created via retrotransposition.
        retrogene,

        // ranscript supported by RNAseq data and not supported by mRNA or EST evidence.
        RNA_Seq_supported_only,

        // ranscript annotated based on mixture of RNA-seq data and EST/mRNA/protein evidence.
        RNA_Seq_supported_partial,

        // ranscript that contains a CDS that has a translation initiation site supported by Ribosomal Profiling data.
        RP_supported_TIS,

        // contains a selenocysteine.
        seleno,

        // a processed pseudogene with one or more introns still present. These are likely formed through the retrotransposition
        semi_processed,

        // ranscript contains at least 1 non-canonical splice junction that is associated with a known or novel genome sequence
        sequence_error,

        // an upstream ATG exists when a downstream ATG is better supported.
        upstream_ATG,

        // a low confidence upstream ATG existing in other coding variant would lead to NMD in this trancript, that uses the high
        upstream_uORF
    }

    public enum TranscriptSupportLevel {
        /** all splice junctions of the transcript are supported by at least one non-suspect mRNA */
        ONE,

        /** the best supporting mRNA is flagged as suspect or the support is from multiple ESTs */
        TWO,

        /** the only support is from a single EST */
        THREE,

        /** the best supporting EST is flagged as suspect */
        FOUR,

        /** no single transcript supports the model structure */
        FIVE,

        /** the transcript was not analyzed */
        NA
    }

    public enum RemapStatus {
        full_contig,
        full_fragment,
        partial,
        deleted,
        no_seq_map,
        gene_conflict,
        gene_size_change,
        automatic_small_ncrna_gene,
        automatic_gene,
        pseudogene
    }

    public enum RemapTargetStatus {
        NEW,
        LOST,
        OVERLAP,
        NONOVERLAP
    }
}
