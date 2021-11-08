import re

from readdatafile.readdatafile import normalize


def readfeaturesfile(featuresfilename: str):
    featurefileheaderlist: list[str] = []
    featurecounter: int
    featurecounter = 0
    featureinfodict: dict = {}
    transcriptgeneassociation: dict = {}
    transcriptconfidence: dict = {}
    transcriptcount: dict = {}
    # 4	module_6	gene_480	transcript_480.6
    # feature id, confidence, count, pathwaymodule, gene, transcript
    # 4483174	1	3	pathway_5	module_67	gene_471	transcript_471.4
    pattern = re.compile("^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+)$")
    with open(featuresfilename) as featurefile:
        # leaves out the line 4 module_6 gene_480 transcript_480.6
        # where confidence and count are unavailable
        featurefilerows = (line.split('\t') for line in featurefile if pattern.match(line))
        for row in featurefilerows:
            if featurecounter != 0:
                loopingcounter: int
                for loopingcounter in range(1, len(featurefileheaderlist)):
                    # row[0] is the feature id
                    # header is _id	confidence	count	Pathway	module	Gene	Transcript
                    # row[0] is transcript information
                    if row[0] not in featureinfodict.keys():
                        featureinfodict[row[0]] = {}
                    if featurefileheaderlist[loopingcounter] not in featureinfodict[row[0]].keys():
                        featureinfodict[row[0]][featurefileheaderlist[loopingcounter]] = row[1:][loopingcounter]
                    # add the transcript - gene association
                    if featurefileheaderlist[loopingcounter] == 'Gene':
                        transcriptgeneassociation[row[0]] = row[1:][loopingcounter]
                    if featurefileheaderlist[loopingcounter] == 'confidence':
                        transcriptconfidence[row[0]] = row[1:][loopingcounter]
                    if featurefileheaderlist[loopingcounter] == 'count':
                        transcriptcount[row[0]] = row[1:][loopingcounter]
            elif featurecounter == 0:
                # header row
                for headerfield in row[1:]:
                    featurefileheaderlist.append(headerfield)
            featurecounter = featurecounter + 1
        maxcount = max(transcriptcount.values())
        transcriptcount = {k: normalize(v, maxcount) for k, v in transcriptcount.items()}
    return featureinfodict, transcriptgeneassociation, transcriptconfidence, transcriptcount