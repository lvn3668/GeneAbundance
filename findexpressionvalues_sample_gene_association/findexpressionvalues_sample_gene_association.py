# Author: Lalitha Viswanathan
# For Inputs
# subjects to yes/no-(treated with drugs)-samples association
# sample to transcripts-per-sample
# transcripts to gene
# mean expression / abundance for each transcript
# returns dict of lists
# key is the gene
# values are list of expression values per transcript for all transcripts associated with that gene

def convertlistoftuplestodict(list_of_tuples: list):
    key_value_dict: dict = {}
    for index, tuple in enumerate(list_of_tuples):
        key_value_dict[tuple[0]] = tuple[1]
        #element_one = tuple[0]
        #element_two = tuple[1]
        #print(element_one, element_two)
    return key_value_dict

def findexpressionvalues_gene_sample_association(subjects_treatedwithdrugs_sample_association: dict,
                                                 sample_to_expressedtranscripts_association: dict,
                                                 transcript_to_gene_association: dict,
                                                 meanexpression_per_sample_for_all_features_or_transcripts_sorted: dict,
                                                 sample_transcript_expressionval_association: dict):


    expressionvalues_subject_drugtreatment_sampleid_transcript_gene: dict = {}
    geneid_expressionvalues: dict = {}
    geneid_to_transcriptexprvals_with_featureids_association: dict = {}

    print(sample_transcript_expressionval_association)
    for subject in subjects_treatedwithdrugs_sample_association.keys():
        if subject not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene.keys():
            expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject] = {}
        feature_or_transcript_treatedwithdrugs_yesno: str
        for feature_or_transcript_treatedwithdrugs_yesno in subjects_treatedwithdrugs_sample_association.get(subject).keys():
            if feature_or_transcript_treatedwithdrugs_yesno not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject].keys():
                expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno] = {}
            for sampleid in subjects_treatedwithdrugs_sample_association.get(subject).get(feature_or_transcript_treatedwithdrugs_yesno):
                transcript_exprvals_dict: dict = {}
                if sampleid in sample_transcript_expressionval_association.keys():
                    transcript_exprvals_dict = convertlistoftuplestodict(sample_transcript_expressionval_association.get(sampleid))
                if sampleid not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno].keys():
                    expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno][sampleid] = {}
                if sampleid in sample_to_expressedtranscripts_association.keys():
                    for transcript in sample_to_expressedtranscripts_association.get(sampleid):
                        if transcript not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno][sampleid].keys():
                            expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno][sampleid][transcript] = {}


                            if transcript_to_gene_association.get(transcript) not in \
                                    expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno][sampleid].keys():

                                expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno][sampleid][transcript][
                                    transcript_to_gene_association.get(transcript)] = {}
                                expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno][sampleid][transcript][
                                    transcript_to_gene_association.get(transcript)] = meanexpression_per_sample_for_all_features_or_transcripts_sorted.get(sampleid)

                                if transcript_to_gene_association.get(transcript) not in geneid_expressionvalues.keys():
                                    geneid_expressionvalues[transcript_to_gene_association.get(transcript)] = []
                                    geneid_to_transcriptexprvals_with_featureids_association[
                                        transcript_to_gene_association.get(transcript)] = None
                                    # Add mean relative abundance to generelativeabundance
                                    # gene expressed via multiple transcripts in multiple samples
                                    # hence list
                                geneid_expressionvalues[transcript_to_gene_association.get(transcript)].append(
                                    meanexpression_per_sample_for_all_features_or_transcripts_sorted.get(sampleid))
                                geneid_to_transcriptexprvals_with_featureids_association[transcript_to_gene_association.get(transcript)] = transcript_exprvals_dict

    return expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues, geneid_to_transcriptexprvals_with_featureids_association