
def normalize(value, divisor=2632):
    return int(value) / int(divisor)



def readdatafile(datafilename: str):
    geneexprvalsperfeaturespersample: dict = {}
    sample_identifier_list: list[str] = []
    featurecounter: int
    featurecounter = 0
    meangeneexpressionpersample: dict = {}
    numberoffeatures: int
    expressedtranscriptspersample : dict = {}
    with open(datafilename) as datafile:
        # split for each sample id, abundence per feature
        rows = (line.split('\t') for line in datafile)
        samplecounter: int
        for row in rows:
            if featurecounter != 0:
                for samplecounter in range(1, len(sample_identifier_list)):
                    if sample_identifier_list[samplecounter] not in meangeneexpressionpersample.keys():
                        meangeneexpressionpersample[sample_identifier_list[samplecounter]] = 0
                    if sample_identifier_list[samplecounter] not in expressedtranscriptspersample.keys():
                        expressedtranscriptspersample[sample_identifier_list[samplecounter]] = []
                    # add abundance measures for each transcript
                    if int(row[1:][samplecounter]) !=0:
                        meangeneexpressionpersample[sample_identifier_list[samplecounter]] = meangeneexpressionpersample.get(
                            sample_identifier_list[samplecounter]) + int(row[1:][samplecounter])
                        # Add expressed transcripts for each sample
                        expressedtranscriptspersample[sample_identifier_list[samplecounter]].append(row[0])
                loopingcounter: int
                for loopingcounter in range(1, len(sample_identifier_list)):
                    # first index is the sample id from the header
                    # inner index is the feature id per row
                    geneexprvalsperfeaturespersample[sample_identifier_list[loopingcounter]][row[0]] = row[1:][loopingcounter]
            elif featurecounter == 0:
                # header row
                # parse to get sample information
                for sample in row[1:]:
                    sample_identifier_list.append(sample)
                    geneexprvalsperfeaturespersample[sample] = {}
            featurecounter = featurecounter + 1

        meangeneexpressionpersample = {k: normalize(v, 2632) for k, v in meangeneexpressionpersample.items()}
        ###########################################################################
        return geneexprvalsperfeaturespersample, meangeneexpressionpersample, expressedtranscriptspersample
