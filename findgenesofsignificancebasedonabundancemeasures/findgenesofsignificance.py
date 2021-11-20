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
    sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength = readfnafile(
        sequencesfilename)

    sample_expressedfeature_or_transcript_expressionvalues_matrix, \
    meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
    expressedtranscriptspersample_sampletotranscriptassoc,sample_transcript_expressionval_association = \
        readsample_features_expressionvalues(
        datafilename)

    sample_expressedfeature_or_transcript_expressionvalues_matrix, samplid_to_subject_association, \
    subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
    sampleid_to_dayssinceexperimentstarted_normalized, sample_treatedwithdrugs_association = readsamplemetadatainformation(
        samplesfilename, sample_expressedfeature_or_transcript_expressionvalues_matrix)

    features_metadata, transcript_or_feature_to_gene_association, transcript_or_feature_to_confidence_normalized_association, \
    transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
    pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
    gene_to_pathway_association, pathway_to_gene_association = readfeaturesmetadatafile(
        featuresfilename)

    sample_to_gene_association, gene_to_sample_association = find_sample_gene_association(
        transcript_or_feature_to_gene_association, expressedtranscriptspersample_sampletotranscriptassoc)

    expressionvalues_subject_drugtreatment_sampleid_transcript_gene: dict
    geneid_expressionvalues: dict
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues, \
    geneid_to_transcriptexprvals_with_featureids_association \
        = findexpressionvalues_gene_sample_association(
        subjects_treatedwithdrugs_to_sample_association,
        expressedtranscriptspersample_sampletotranscriptassoc,
        transcript_or_feature_to_gene_association,
        meanexpression_per_sample_for_all_features_or_transcripts_sorted,sample_transcript_expressionval_association)

    return sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength, \
           sample_expressedfeature_or_transcript_expressionvalues_matrix, \
           meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
           expressedtranscriptspersample_sampletotranscriptassoc, sample_expressedfeature_or_transcript_expressionvalues_matrix, \
           samplid_to_subject_association, \
           subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
           sampleid_to_dayssinceexperimentstarted_normalized, features_metadata, transcript_or_feature_to_gene_association, \
           transcript_or_feature_to_confidence_normalized_association, \
           transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
           pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
           gene_to_pathway_association, pathway_to_gene_association, sample_to_gene_association, gene_to_sample_association, \
           expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues, \
           sample_treatedwithdrugs_association, sample_transcript_expressionval_association, geneid_to_transcriptexprvals_with_featureids_association
