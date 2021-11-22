# Author: Lalitha Viswanathan
# Package to read metadata (samples, features and expression data)
from findexpressionvalues_sample_gene_association.findexpressionvalues_sample_gene_association import \
    findexpressionvalues_gene_sample_association
from readdatafile.readdatafile import readsample_features_expressionvalues
from readfeaturefile.readfeaturesfiles import readfeaturesmetadatafile
from readsamplemetadatafile.readsamplemetadatainformation import readsamplemetadatainformation
from readsequencefiles.readsequencesfna import readfnafile
from sample_to_gene_association.sample_gene_association import find_sample_gene_association


def readmetadataaboutsamplesfeaturesandexpressionvalues(datafilename: str, featuresfilename: str, samplesfilename: str,
                                                        sequencesfilename: str):
    ### read sequence information
    sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength = readfnafile(
        sequencesfilename)
#######################################################################################################

    # read sample metadata and return
    # expression values
    # mean expr values per sample
    # expressed transcripts per sample
    # samples containing a given transcript (per transcript)
    sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals, \
    meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
    expressedtranscriptspersample_sampletotranscriptassoc_noexprvals, persample_expressed_transcript_to_expressionval_association, \
    pertranscript_sampleidentifier_expressionvals_across_samples = \
        readsample_features_expressionvalues(
            datafilename)

#######################################################################################################

    # return matrix populated with samples and transcripts
    # Sample to subject association
    # subjects-treated with drugs to sample association
    # samples with no transcripts and hence no expression measures
    # samples to days in experimental pool
    # samples to treated with drugs association
    sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals, samplid_to_subject_association, \
    subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
    sampleid_to_dayssinceexperimentstarted_normalized, sample_treatedwithdrugs_association = readsamplemetadatainformation(
        samplesfilename, sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals)

#######################################################################################################

    # features metadat file
    # transcript to gene association
    # transcript confidence normalized
    # transcript count normalized
    # pathway ranking based on number of modules per pathway and number of genes in each module
    # gene to pathway association
    # pathway to gene association
    features_metadata, transcript_or_feature_to_gene_association, transcript_or_feature_to_confidence_normalized_association, \
    transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
    pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
    gene_to_pathway_association, pathway_to_gene_association = readfeaturesmetadatafile(
        featuresfilename)
#######################################################################################################
    # sample to gene association
    # gene to sample association
    sample_to_gene_association, gene_to_sample_association = find_sample_gene_association(
        transcript_or_feature_to_gene_association, expressedtranscriptspersample_sampletotranscriptassoc_noexprvals)

#######################################################################################################
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene: dict
    geneid_expressionvalues: dict
    geneid_to_all_transcriptexprvals_acrosssamples_association = \
        findexpressionvalues_gene_sample_association(
            subjects_treatedwithdrugs_to_sample_association,
            transcript_or_feature_to_gene_association,
            pertranscript_sampleidentifier_expressionvals_across_samples,
            gene_to_transcript_or_feature_association, transcript_or_feature_to_confidence_normalized_association, \
            transcript_or_feature_to_count_normalized_association, sampleid_to_dayssinceexperimentstarted_normalized,
            sampleid_to_length_association_normalizedoverrmaxlength)

    return sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength, \
sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals, \
    meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
    expressedtranscriptspersample_sampletotranscriptassoc_noexprvals, persample_expressed_transcript_to_expressionval_association, \
    pertranscript_sampleidentifier_expressionvals_across_samples, samplid_to_subject_association, \
    subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
    sampleid_to_dayssinceexperimentstarted_normalized, sample_treatedwithdrugs_association, features_metadata, transcript_or_feature_to_gene_association, transcript_or_feature_to_confidence_normalized_association, \
    transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
    pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
    gene_to_pathway_association, pathway_to_gene_association, sample_to_gene_association, gene_to_sample_association, geneid_to_all_transcriptexprvals_acrosssamples_association


