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

# Function reads expression values matrix and returns
# per sample:
# expressed transcript list
# expressed transcript and expression value
# per transcript, sample id and expression value (across samples)
# mean of expression measures across transcripts (per sample) ;
# mean calculated as sum of non-zero expression measures for each transcript per sample
# divide sum by number of expressed transcripts per sample
def readsample_features_expressionvalues(expressionvaluesfilename: str):
    sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals: dict = {}
    sample_identifier_list: list[str] = []
    featurecounter: int
    featurecounter = 0
    meanexpression_per_sample_for_all_expressed_features_or_transcripts: dict = {}
    numberoffeatures: int
    expressedtranscriptspersample_sampletotranscriptassoc_noexprvals: dict = {}
    persample_expressed_transcript_to_expressionval_association: dict = {}
    pertranscript_sampleidentifier_expressionvals_across_samples: dict = {}
    with open(expressionvaluesfilename) as datafile:
        # split for each sample id, abundence per feature
        rows = (line.split('\t') for line in datafile)
        samplecounter: int
        featureortranscriptrow: list[str]
        for featureortranscriptrow in rows:
            if featurecounter != 0:
                if featureortranscriptrow[0] not in pertranscript_sampleidentifier_expressionvals_across_samples.keys():
                    pertranscript_sampleidentifier_expressionvals_across_samples[featureortranscriptrow[0].strip()] = []
                for samplecounter in range(1, len(sample_identifier_list)):
                    if sample_identifier_list[
                        samplecounter] not in meanexpression_per_sample_for_all_expressed_features_or_transcripts.keys():
                        meanexpression_per_sample_for_all_expressed_features_or_transcripts[
                            sample_identifier_list[samplecounter]] = 0
                    if sample_identifier_list[
                        samplecounter] not in expressedtranscriptspersample_sampletotranscriptassoc_noexprvals.keys():
                        expressedtranscriptspersample_sampletotranscriptassoc_noexprvals[
                            sample_identifier_list[samplecounter]] = []
                        persample_expressed_transcript_to_expressionval_association[
                            sample_identifier_list[samplecounter]] = []

                    # add abundance measures for each transcript
                    if int(featureortranscriptrow[1:][samplecounter]) != 0:
                        meanexpression_per_sample_for_all_expressed_features_or_transcripts[
                            sample_identifier_list[
                                samplecounter]] = meanexpression_per_sample_for_all_expressed_features_or_transcripts.get(
                            sample_identifier_list[samplecounter]) + int(featureortranscriptrow[1:][samplecounter])

                        # Add expressed transcripts/ features id for each sample
                        expressedtranscriptspersample_sampletotranscriptassoc_noexprvals[
                            sample_identifier_list[samplecounter]].append(featureortranscriptrow[0])
                        persample_expressed_transcript_to_expressionval_association[
                            sample_identifier_list[samplecounter]].append(
                            (featureortranscriptrow[0], int(featureortranscriptrow[1:][samplecounter])))

                loopingcounter: int
                for loopingcounter in range(1, len(sample_identifier_list)):
                    # first index is the sample id from the header
                    # inner index is the feature id per featureortranscriptrow
                    sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals[
                        sample_identifier_list[loopingcounter]][featureortranscriptrow[0]] = featureortranscriptrow[1:][
                        loopingcounter]
                    if int(featureortranscriptrow[1:][loopingcounter]) != 0:
                        pertranscript_sampleidentifier_expressionvals_across_samples[
                            featureortranscriptrow[0].strip()].append((sample_identifier_list[loopingcounter].strip(),
                                                                       int(featureortranscriptrow[1:][loopingcounter])))
            elif featurecounter == 0:
                # header featureortranscriptrow
                # parse to get sample information
                for sample in featureortranscriptrow[1:]:
                    sample_identifier_list.append(sample)
                    sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals[sample] = {}
            featurecounter = featurecounter + 1

        # mean is normalized over all features or transcripts per sample
        meanexpression_per_sample_for_all_expressed_features_or_transcripts: dict = {
            k: normalize(v, len(expressedtranscriptspersample_sampletotranscriptassoc_noexprvals.get(k))) for k, v in
            meanexpression_per_sample_for_all_expressed_features_or_transcripts.items()}

        meanexpression_per_sample_for_all_features_or_transcripts_sorted = dict(
            sorted(meanexpression_per_sample_for_all_expressed_features_or_transcripts.items(),
                   key=lambda item: item[1],
                   reverse=True))

        ###########################################################################
        persample_expressed_transcript_to_expressionval_association[sample_identifier_list[samplecounter]].append(
            (featureortranscriptrow[0], int(featureortranscriptrow[1:][samplecounter])))

        ###########################################################################
        return sample_expressedfeature_or_transcript_expressionvalues_matrix_allvals, \
               meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
               expressedtranscriptspersample_sampletotranscriptassoc_noexprvals, persample_expressed_transcript_to_expressionval_association, \
               pertranscript_sampleidentifier_expressionvals_across_samples
