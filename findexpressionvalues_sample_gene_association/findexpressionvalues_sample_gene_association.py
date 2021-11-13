def findexpressionvalues_gene_sample_association(subjects_treatedwithdrugs_sample_association: dict, sample_to_expressedtranscripts_association: dict,
                                                 transcript_to_gene_association: dict, meanexpression_per_sample_for_all_features_or_transcripts_sorted: dict):


    expressionvalues_subject_drugtreatment_sampleid_transcript_gene: dict = {}
    geneid_expressionvalues: dict = {}
    for subject in subjects_treatedwithdrugs_sample_association.keys():
        if subject not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene.keys():
            expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject] = {}
        feature_or_transcript_treatedwithdrugs_yesno: str
        for feature_or_transcript_treatedwithdrugs_yesno in subjects_treatedwithdrugs_sample_association.get(subject).keys():
            if feature_or_transcript_treatedwithdrugs_yesno not in expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject].keys():
                expressionvalues_subject_drugtreatment_sampleid_transcript_gene[subject][feature_or_transcript_treatedwithdrugs_yesno] = {}
            for sampleid in subjects_treatedwithdrugs_sample_association.get(subject).get(feature_or_transcript_treatedwithdrugs_yesno):
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
                                    # Add mean relative abundance to generelativeabundance
                                    # gene expressed via multiple transcripts in multiple samples
                                    # hence list
                                geneid_expressionvalues[transcript_to_gene_association.get(transcript)].append(
                                    meanexpression_per_sample_for_all_features_or_transcripts_sorted.get(sampleid))
    return expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues