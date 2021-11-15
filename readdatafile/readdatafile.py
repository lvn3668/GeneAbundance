# Author Lalitha Viswanathan
# Gene Abundance
# Read expression values per sample (total of 35 samples) and 2000-odd features / transcripts

def normalize(value, divisor):
    if float(divisor) != 0.0:
        return float(value) / float(divisor)
    else:
        return float(divisor)


# read data file of expression values per sample per transcript
# not all transcripts are expressed in each sample
# returns the matrix
# normalizes the expression values per transcript per sample
# over number of transcripts per sample
#
# in addition returns number of transcripts per sample (dict of lists)
def readsampletofeatureexpressionvalues(expressionvaluesfilename: str):
    sample_expressedfeature_or_transcript_expressionvalues_matrix: dict = {}
    sample_identifier_list: list[str] = []
    featurecounter: int
    featurecounter = 0
    meanexpression_per_sample_for_all_features_or_transcripts: dict = {}
    numberoffeatures: int
    expressedtranscriptspersample_sampletotranscriptassoc: dict = {}
    with open(expressionvaluesfilename) as datafile:
        # split for each sample id, abundence per feature
        rows = (line.split('\t') for line in datafile)
        samplecounter: int
        featureortranscriptrow: list[str]
        for featureortranscriptrow in rows:
            if featurecounter != 0:
                for samplecounter in range(1, len(sample_identifier_list)):
                    if sample_identifier_list[
                        samplecounter] not in meanexpression_per_sample_for_all_features_or_transcripts.keys():
                        meanexpression_per_sample_for_all_features_or_transcripts[
                            sample_identifier_list[samplecounter]] = 0
                    if sample_identifier_list[
                        samplecounter] not in expressedtranscriptspersample_sampletotranscriptassoc.keys():
                        expressedtranscriptspersample_sampletotranscriptassoc[
                            sample_identifier_list[samplecounter]] = []

                    # add abundance measures for each transcript
                    if int(featureortranscriptrow[1:][samplecounter]) != 0:
                        meanexpression_per_sample_for_all_features_or_transcripts[
                            sample_identifier_list[
                                samplecounter]] = meanexpression_per_sample_for_all_features_or_transcripts.get(
                            sample_identifier_list[samplecounter]) + int(featureortranscriptrow[1:][samplecounter])

                        # Add expressed transcripts/ features id for each sample
                        expressedtranscriptspersample_sampletotranscriptassoc[
                            sample_identifier_list[samplecounter]].append(featureortranscriptrow[0])
                loopingcounter: int
                for loopingcounter in range(1, len(sample_identifier_list)):
                    # first index is the sample id from the header
                    # inner index is the feature id per featureortranscriptrow
                    sample_expressedfeature_or_transcript_expressionvalues_matrix[
                        sample_identifier_list[loopingcounter]][featureortranscriptrow[0]] = featureortranscriptrow[1:][
                        loopingcounter]
            elif featurecounter == 0:
                # header featureortranscriptrow
                # parse to get sample information
                for sample in featureortranscriptrow[1:]:
                    sample_identifier_list.append(sample)
                    sample_expressedfeature_or_transcript_expressionvalues_matrix[sample] = {}
            featurecounter = featurecounter + 1

        # normalize over number of expressed transcripts / features per sample
        meanexpression_per_sample_for_all_features_or_transcripts: dict = {
            k: normalize(v, len(expressedtranscriptspersample_sampletotranscriptassoc.get(k))) for k, v in
            meanexpression_per_sample_for_all_features_or_transcripts.items()}

        meanexpression_per_sample_for_all_features_or_transcripts_sorted = dict(
            sorted(meanexpression_per_sample_for_all_features_or_transcripts.items(),
                   key=lambda item: item[1],
                   reverse=True))

        ###########################################################################
        return sample_expressedfeature_or_transcript_expressionvalues_matrix, \
               meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
               expressedtranscriptspersample_sampletotranscriptassoc
