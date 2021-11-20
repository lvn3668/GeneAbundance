# Author Lalitha Viswanathan
# Function to normalize expression values per transcript
# over confidence measure, count measure
# sequence length normalized
# pathway rank

def normalizeexpressionmeasures(gene: str, geneid_expressionvalues: dict,
                                genesrankedbypathwayexpression: dict,
                                gene_to_transcript_or_feature_association: dict,
                                transcript_or_feature_to_confidence_association: dict,
                                transcript_or_feature_to_count_association: dict,
                                gene_to_sample_association: dict,
                                sampleid_to_dayssinceexperimentstarted_normalized: dict,
                                sampleid_to_length_association_normalizedoverrmaxlength: dict,
                                genesmissingpathwayassociations: list) -> list:
    """

    :type genesrankedbypathwayexpression: object
    :param gene:
    :param geneid_expressionvalues:
    :param genesrankedbypathwayexpression:
    :param gene_to_transcript_or_feature_association:
    :param transcript_or_feature_to_confidence_association:
    :param transcript_or_feature_to_count_association:
    :param gene_to_sample_association:
    :param sampleid_to_dayssinceexperimentstarted_normalized:
    :param sampleid_to_length_association_normalizedoverrmaxlength:
    :param genesmissingpathwayassociations:
    :return:
    """
    listofabundancevals: list[str] = []

    for expressionvaluesnormalizedoverexpressedfeatures_or_transcripts in geneid_expressionvalues.get(gene):
        if genesrankedbypathwayexpression.get(gene) != 0:
            expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts * genesrankedbypathwayexpression.get(
                gene)
        else:
            genesmissingpathwayassociations.append(gene)
        # key is gene
        # value is normalized expression values
        sum_transcript_count_measure: float
        sum_transcript_count_measure = 0
        sum_transcript_confidence_measure: float
        sum_transcript_confidence_measure = 0
        for transcript in gene_to_transcript_or_feature_association.get(gene):


            # sum all the transcript confidence measures
            sum_transcript_confidence_measure = 0
            if transcript in transcript_or_feature_to_confidence_association.keys():
                sum_transcript_confidence_measure = sum_transcript_confidence_measure + float(
                    transcript_or_feature_to_confidence_association.get(transcript))

            # then normalize by count
            sum_transcript_count_measure = 0
            if transcript in transcript_or_feature_to_count_association.keys():
                sum_transcript_count_measure = sum_transcript_count_measure + float(
                    transcript_or_feature_to_count_association.get(transcript))

        # normalize by confidence and count (positive correlation)
        expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts * sum_transcript_confidence_measure * sum_transcript_count_measure

        # abundance values normalized by count and confidence level
        # Use days since experiment started for treated-with-drug NO
        # If it is NOT treated with drugs and days since experiment started is greater than 0
        # Inverse proportion
        # Longer the sample in the experiment, and not treated with drug, over expressed it is,
        # Normalize days since experiment started (0/ (max), 10 /(max) , etc.)
        # Abundance is unaltered for treated-with-drug yes cases (????)
        # For untreated cases, account for over expression / under expression by DIVIDING the abundance measure
        # by normalized days-since-experiment-started
        # This will underexpress those transcripts belonging to UNTREATED samples,
        # that have been the longest in the experiment pool
        sumdayssinceexperimentstarted: float
        sumdayssinceexperimentstarted = 0
        for sample in gene_to_sample_association.get(gene):
            sumdayssinceexperimentstarted = 0
            if sample in sampleid_to_dayssinceexperimentstarted_normalized.keys() and \
                    sampleid_to_dayssinceexperimentstarted_normalized.get(sample) != 0:
                sumdayssinceexperimentstarted = sumdayssinceexperimentstarted + sampleid_to_dayssinceexperimentstarted_normalized.get(
                    sample)

        if sumdayssinceexperimentstarted != 0:
            expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts / sumdayssinceexperimentstarted

        # gene is part of multiple samples
        # for each sample, that the gene is expressed in
        # get sequence length (normalized over longest sequence length)
        # aggregate normalized sequence lengths
        # divide by the sum
        # negative correlation (longer the sequence, longer the PCR fragments and higher expression values)
        # correction is by dividing the sum of lengths
        sumnormalizedsequencelengths: float
        sumnormalizedsequencelengths = 0
        for sample in gene_to_sample_association.get(gene):
            sumnormalizedsequencelengths = 0
            if sample in sampleid_to_length_association_normalizedoverrmaxlength.keys() and \
                    sampleid_to_length_association_normalizedoverrmaxlength.get(sample) != 0:
                sumnormalizedsequencelengths = sumnormalizedsequencelengths + sampleid_to_length_association_normalizedoverrmaxlength.get(
                    sample)

        if sumnormalizedsequencelengths != 0:
            expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts / sumnormalizedsequencelengths
        listofabundancevals.append(expressionvaluesnormalizedoverexpressedfeatures_or_transcripts)

    # return dict of transcript id and abundance vals
    return listofabundancevals