# Author Lalitha Viswanathan
# Function to normalize expression values per transcript
# over confidence measure, count measure
# sequence length normalized
# pathway rank
from readdatafile.readdatafile import normalize


def normalizeexpressionmeasures(gene:str,  transcriptandexprmeasure_unnormalized: dict,
                                  sampleid_to_dayssinceexperimentstarted_normalized:dict,
                                  sampleid_to_length_association_normalizedoverrmaxlength:dict,
                                  genesrankedbypathwayexpression: dict
                                  ) -> (list, dict):
    """

    :type genesrankedbypathwayexpression: object
    :param gene:
    :param geneid_expressionvalues:
    :param genesrankedbypathwayexpression:
    :param gene_to_transcript_or_feature_association:
    :param transcript_or_feature_to_confidence_normalized_association:
    :param transcript_or_feature_to_count_normalized_association:
    :param gene_to_sample_association:
    :param sampleid_to_dayssinceexperimentstarted_normalized:
    :param sampleid_to_length_association_normalizedoverrmaxlength:
    :param genesmissingpathwayassociations:
    :return:
    """

    per_transcript_expression_measure_tobenormalized_across_samples: dict = {
        k: normalize(v, max(transcriptandexprmeasure_unnormalized.values())) for k, v in
        transcriptandexprmeasure_unnormalized.items()}

    print(per_transcript_expression_measure_tobenormalized_across_samples)
    per_transcript_expression_measure_tobenormalized_across_samples = dict(
        sorted(per_transcript_expression_measure_tobenormalized_across_samples.items(),
               key=lambda item: item[1],
               reverse=True))
    # multiply by pathway rank
    #per_transcript_expression_measure_tobenormalized_across_samples = {key: value * genesrankedbypathwayexpression.get(
    #            gene)  & print(key, value) for key, value in per_transcript_expression_measure_tobenormalized_across_samples.items()
    #                                                                   if int(genesrankedbypathwayexpression.get(
    #            gene)) != 0}
    #print("After pathway rank")
    #print(per_transcript_expression_measure_tobenormalized_across_samples)

    # divide by days since in experimental pool
    per_transcript_expression_measure_tobenormalized_across_samples = \
        {key: value / (sampleid_to_dayssinceexperimentstarted_normalized.get(
                    key))  if key in sampleid_to_dayssinceexperimentstarted_normalized.keys()
                        and int(sampleid_to_dayssinceexperimentstarted_normalized.get(key)) !=0 else value for key, value in per_transcript_expression_measure_tobenormalized_across_samples.items()  }

    print("After days in exp pool")
    print(per_transcript_expression_measure_tobenormalized_across_samples)

    # divide by sequence length
    per_transcript_expression_measure_tobenormalized_across_samples = {
        key: value / sampleid_to_length_association_normalizedoverrmaxlength.get(
            key) if key in sampleid_to_length_association_normalizedoverrmaxlength.keys() else value for key, value in per_transcript_expression_measure_tobenormalized_across_samples.items()
    }

    ################################################################################

    # abundance values normalized by count and confidence level
    # Use days since experiment started for treated-with-drug NO
    # If it is NOT treated with drugs and days since experiment started is greater than 0
    # Inverse proportion
    #  Longer the sample in the experiment, and not treated with drug, over expressed it is,
    # Normalize days since experiment started (0/ (max), 10 /(max) , etc.)
    # Abundance is unaltered for treated-with-drug yes cases (????)
    # For untreated cases, account for over expression / under expression by DIVIDING the abundance measure
    # by normalized days-since-experiment-started
    # This will underexpress those transcripts belonging to UNTREATED samples,
    # that have been the longest in the experiment pool
        # gene is part of multiple samples
        # for each sample, that the gene is expressed in
        # get sequence length (normalized over longest sequence length)
        # aggregate normalized sequence lengths
        # divide by the sum
        # negative correlation (longer the sequence, longer the PCR fragments and higher expression values)
        # correction is by dividing the sum of lengths

    # return dict of transcript id and abundance vals
    print(per_transcript_expression_measure_tobenormalized_across_samples)
    #exit(1)
    return per_transcript_expression_measure_tobenormalized_across_samples