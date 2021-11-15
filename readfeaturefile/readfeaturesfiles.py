# Author Lalitha Viswanathan
# Program to read features meta data such as counts, confidence measure, pathway, module within pathway, gene

import re
from typing import TextIO
from pathwayranking.pathwayrank import pathwayrank
from datanormalization.normalizeovermaxvalue import normalizeovermaximum
# read features metadata file
# assumption of records containing all the fields  (namely id, count, confidence, pathway, gene, transcripts)
# other records are ignored
def readfeaturesmetadatafile(featuresfilename: str) -> tuple[dict, dict, dict, dict, dict, dict, dict, dict]:
    featurefileheaderlist: list[str] = []
    featurecounter: int
    featurecounter = 0
    featureinfodict: dict = {}
    transcript_or_feature_to_gene_association: dict = {}
    transcript_or_feature_to_confidence_normalized_association: dict = {}
    transcript_or_feature_to_count_normalized_association: dict = {}
    gene_to_transcript_or_feature_association: dict = {}
    pathway_to_module_nummodulesperpathway_association: dict = {}
    gene_to_pathway_association: dict = {}
    pathway_to_gene_association: dict = {}
    # 4	module_6	gene_480	transcript_480.6
    # feature id, confidence, count, pathwaymodule, gene, transcript
    # 4483174	1	3	pathway_5	module_67	gene_471	transcript_471.4
    pattern = re.compile("^(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+)$")
    featurefile: TextIO
    with open(featuresfilename) as featurefile:
        # leaves out the line 4 module_6 gene_480 transcript_480.6
        # where confidence and count are unavailable
        featurefilerows = (line.split('\t') for line in featurefile if pattern.match(line))
        for row in featurefilerows:
            if featurecounter != 0:
                counter: int
                for counter in range(0, len(featurefileheaderlist)):
                    # print("Feature counter ", featurecounter, " row[0] ", row[0], " featurefileheaderlist[counter] ",
                    #      featurefileheaderlist[counter])
                    # Populate the featureinfodict with feature-id and metadatafield (count / confidence / pathway / module etc.)
                    # as inner key
                    # Thus dict entries are of the form [id][count] = [value], [id][confidence] = value etc.
                    ######################################################################
                    if row[0] not in featureinfodict.keys():
                        featureinfodict[row[0]] = {}

                    if featurefileheaderlist[counter] not in featureinfodict[row[0]].keys():
                        featureinfodict[row[0]][featurefileheaderlist[counter]] = row[1:][counter]
                    #######################################################################
                    # add the transcript - gene association
                    # Key is feature id
                    # Value is associated gene
                    # Transcript Id as defined in features file is not captured
                    # as abundance measures refer to feature id
                    # transcript - feature not captured ; hence used interchangeably
                    if featurefileheaderlist[counter] == 'Gene':
                        transcript_or_feature_to_gene_association[row[0]] = row[1:][counter]
                        #######################################################################
                        # Add all transcripts (noted as feature id) associated with a given gene
                        if row[1:][counter] not in gene_to_transcript_or_feature_association.keys():
                            gene_to_transcript_or_feature_association[row[1:][counter]] = []
                        # Append transcript / feature information
                        gene_to_transcript_or_feature_association[row[1:][counter]].append(row[0])

                    #######################################################################
                    # Add confidence measure associated with given feature
                    if featurefileheaderlist[counter] == 'confidence':
                        transcript_or_feature_to_confidence_normalized_association[row[0]] = row[1:][counter]
                    #######################################################################
                    # Add counts associated with given feature
                    if featurefileheaderlist[counter] == 'count':
                        transcript_or_feature_to_count_normalized_association[row[0]] = row[1:][counter]
                    #######################################################################
                    # Add pathway ranking measure
                    # Add association bewteen pathway and num modules ;
                    # pathway and gene (1: many association)
                    # gene and pathway (1: many association)
                    if featurefileheaderlist[counter] == 'Pathway':
                        # check if pathway is in the dict
                        # create pathway module and pathway gene association
                        if row[1:][counter] not in pathway_to_module_nummodulesperpathway_association.keys():
                            pathway_to_module_nummodulesperpathway_association[row[1:][counter]] = []
                        if row[1:][counter] not in pathway_to_gene_association.keys():
                            pathway_to_gene_association[row[1:][counter]] = []
                        pathway_to_module_nummodulesperpathway_association[row[1:][counter]].append(
                            row[1:][counter + 1])
                        pathway_to_gene_association[row[1:][counter]].append(row[1:][counter + 2])
                        # check if gene is in the dict
                        if row[1:][counter + 2] not in gene_to_pathway_association.keys():
                            gene_to_pathway_association[row[1:][counter + 2]] = []
                        gene_to_pathway_association[row[1:][counter + 2]].append(row[1:][counter])
            elif featurecounter == 0:
                # header row
                headerfield: str
                for headerfield in row[1:]:
                    featurefileheaderlist.append(headerfield)
                print("Feature row header row ", featurefileheaderlist)
            featurecounter = featurecounter + 1
        featurefile.close()
        ###################################################################################################
        # Normalize the featureid to count association
        transcript_or_feature_to_count_normalized_association = normalizeovermaximum(transcript_or_feature_to_count_normalized_association)
        ###################################################################################################
        # Normalize the featureid to confidence association
        transcript_or_feature_to_confidence_normalized_association = normalizeovermaximum(transcript_or_feature_to_confidence_normalized_association)
        ###################################################################################################
        # print(json.dumps(gene_to_transcript_or_feature_association, sort_keys=True, indent=2))
        # print(json.dumps(pathway_to_module_nummodulesperpathway_association, sort_keys=True, indent=2))
        # print(json.dumps(gene_to_pathway_association, sort_keys=True, indent=2))
        # Normalize pathway gene association over number of genes per pathway
        pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes = pathwayrank(pathway_to_gene_association,
                                                                                    pathway_to_module_nummodulesperpathway_association)
    # Values returned
    # Feature metadata dict
    # feature to gene association
    # feature to confidence association
    # feature to count normalized association
    # gene to transcript normalized
    # pathway rank normalized by weighted gene-numberOfModules correlation
    # gene to pathway association
    # pathway to gene association
    return featureinfodict, transcript_or_feature_to_gene_association, \
           transcript_or_feature_to_confidence_normalized_association, \
           transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
           pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
           gene_to_pathway_association, pathway_to_gene_association
