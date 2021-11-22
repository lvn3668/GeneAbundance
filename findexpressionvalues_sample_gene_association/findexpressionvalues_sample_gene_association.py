# Author: Lalitha Viswanathan
# For Inputs
# subjects to yes/no-(treated with drugs)-samples association
# sample to transcripts-per-sample
# transcripts to gene
# mean expression / abundance for each transcript
# returns dict of lists
# key is the gene
# values are list of expression values per transcript for all transcripts associated with that gene

def normalizebycountandconfidence(persampleexprvalsforeachtranscript: dict,
                                  transcript: str,
                                  transcript_or_feature_to_confidence_normalized_association: dict,
                                  transcript_or_feature_to_count_normalized_association: dict,
                                  sampleid_to_dayssinceexperimentstarted_normalized: dict,
                                  sampleid_to_length_association_normalizedoverrmaxlength: dict
                                  ):
    # normalize by confidence and count
    persampleexprvalsforeachtranscript = {
        key: value * transcript_or_feature_to_confidence_normalized_association.get(transcript) *
             transcript_or_feature_to_count_normalized_association.get(transcript)
        for key, value in persampleexprvalsforeachtranscript.items()}
    for k, v in persampleexprvalsforeachtranscript.items():
        if k in sampleid_to_length_association_normalizedoverrmaxlength.keys() and \
                sampleid_to_length_association_normalizedoverrmaxlength.get(k) != 0:
            persampleexprvalsforeachtranscript[k] = persampleexprvalsforeachtranscript.get(
                k) / sampleid_to_length_association_normalizedoverrmaxlength.get(k)

        if k in sampleid_to_dayssinceexperimentstarted_normalized.keys() and \
                sampleid_to_dayssinceexperimentstarted_normalized.get(k) != 0:
            persampleexprvalsforeachtranscript[k] = persampleexprvalsforeachtranscript.get(
                k) / sampleid_to_dayssinceexperimentstarted_normalized.get(k)

    return persampleexprvalsforeachtranscript


def convertlistoftuplestodict(list_of_tuples: list):
    key_value_dict: dict = {}
    for index, tuple in enumerate(list_of_tuples):
        key_value_dict[tuple[0]] = tuple[1]
    return key_value_dict


def findexpressionvalues_gene_sample_association(
        subjects_treatedwithdrugs_to_sample_association: dict,
            transcript_or_feature_to_gene_association: dict,
            pertranscript_sampleidentifier_expressionvals_across_samples: dict,
            gene_to_transcript_or_feature_association: dict,
                                                 transcript_or_feature_to_confidence_normalized_association: dict,
            transcript_or_feature_to_count_normalized_association: dict,
                                                 sampleid_to_dayssinceexperimentstarted_normalized: dict,
            sampleid_to_length_association_normalizedoverrmaxlength: dict):
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene: dict = {}
    geneid_expressionvalues_foralltrancripts_persample: dict = {}
    geneid_to_all_transcriptexprvals_acrosssamples_association: dict = {}

    # for each gene convert dicts to counters and sum them
    from collections import Counter

    for gene in gene_to_transcript_or_feature_association.keys():
        geneid_to_all_transcriptexprvals_acrosssamples_association[gene] = []
        for transcript in gene_to_transcript_or_feature_association.get(gene):
            # transcript to gene is 1:1 association passes sanity check
            # can be present in multiple samples
            # A gene can have multiple trancripts
            # so add transcript measure across samples for each transcript
            geneid_to_all_transcriptexprvals_acrosssamples_association[gene].append(
                normalizebycountandconfidence(
                    convertlistoftuplestodict(
                        pertranscript_sampleidentifier_expressionvals_across_samples.get(transcript.strip())),
                    transcript, transcript_or_feature_to_confidence_normalized_association,
                    transcript_or_feature_to_count_normalized_association,
                    sampleid_to_dayssinceexperimentstarted_normalized,
                    sampleid_to_length_association_normalizedoverrmaxlength
                ))
        z: dict = {}
        for d in geneid_to_all_transcriptexprvals_acrosssamples_association.get(gene):
            z = dict(Counter(d) + Counter(z))
        geneid_to_all_transcriptexprvals_acrosssamples_association[gene] = dict(z)
        print("Gene ", gene, " xscript expr vals aggregated ",
              geneid_to_all_transcriptexprvals_acrosssamples_association[gene])

    # for subject in subjects_treatedwithdrugs_to_sample_association.keys():
    #     if subject not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene.keys():
    #         expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject] = {}
    #     feature_or_transcript_treatedwithdrugs_yesno: str
    #     for feature_or_transcript_treatedwithdrugs_yesno in subjects_treatedwithdrugs_to_sample_association.get(
    #             subject).keys():
    #         if feature_or_transcript_treatedwithdrugs_yesno not in \
    #                 expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject].keys():
    #             expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                 feature_or_transcript_treatedwithdrugs_yesno] = {}
    #         for sampleid in subjects_treatedwithdrugs_to_sample_association.get(subject).get(
    #                 feature_or_transcript_treatedwithdrugs_yesno):
    #             alltranscripts_insample_to_corresponding_exprvals_dict: dict = {}
    #             if sampleid in subjects_treatedwithdrugs_to_sample_association.keys():
    #                 # get list of transcripts and expr vals (unnormalized) per sample
    #                 alltranscripts_insample_to_corresponding_exprvals_dict = convertlistoftuplestodict(
    #                     pertranscript_sampleidentifier_expressionvals_across_samples.get(sampleid))
    #             if sampleid not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                 feature_or_transcript_treatedwithdrugs_yesno].keys():
    #                 expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                     feature_or_transcript_treatedwithdrugs_yesno][sampleid] = {}
    #             if sampleid in sample_to_expressedtranscripts_association.keys():
    #                 for transcript in sample_to_expressedtranscripts_association.get(sampleid):
    #
    #                     if transcript not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                         feature_or_transcript_treatedwithdrugs_yesno][sampleid].keys():
    #                         expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                             feature_or_transcript_treatedwithdrugs_yesno][sampleid][transcript] = {}
    #
    #                         if transcript_to_gene_association.get(transcript) not in \
    #                                 expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                                     feature_or_transcript_treatedwithdrugs_yesno][sampleid].keys():
    #
    #                             expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                                 feature_or_transcript_treatedwithdrugs_yesno][sampleid][transcript][
    #                                 transcript_to_gene_association.get(transcript)] = {}
    #                             expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][
    #                                 feature_or_transcript_treatedwithdrugs_yesno][sampleid][transcript][
    #                                 transcript_to_gene_association.get(
    #                                     transcript)] = meanexpression_per_sample_for_all_features_or_transcripts_sorted.get(
    #                                 sampleid)
    #
    #                             if transcript_to_gene_association.get(
    #                                     transcript) not in geneid_expressionvalues_foralltrancripts_persample.keys():
    #                                 geneid_expressionvalues_foralltrancripts_persample[
    #                                     transcript_to_gene_association.get(transcript)] = []
    #
    #                                 # Add mean relative abundance to generelativeabundance
    #                                 # gene expressed via multiple transcripts in multiple samples
    #                                 # hence list
    #                                 # mean expresion values per sample that the gene's transcripts are expressed in
    #                             geneid_expressionvalues_foralltrancripts_persample[
    #                                 transcript_to_gene_association.get(transcript)].append(
    #                                 meanexpression_per_sample_for_all_features_or_transcripts_sorted.get(sampleid))

    # sub - sample - trancript - gene
    # at this stage, the dict geneid_to_all_transcriptexprvals_acrosssamples_association is
    # mapping between gene and transcript-expr vals (dict)
    # starting with subject, the associated sample
    # from sample, associated transcripts
    # from transcript , associated gene (1:1)
    # for each sample, subdict of transcripts - expr vals
    # store gene - subdict association (multiple)
    return geneid_to_all_transcriptexprvals_acrosssamples_association
