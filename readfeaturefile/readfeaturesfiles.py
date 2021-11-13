# Author Lalitha Viswanathan
# Program to read features meta data such as counts, confidence intervals, pathway, module within pathway, gene

import re
from typing import TextIO

import numpy
import scipy
import scipy.stats
from readdatafile.readdatafile import normalize
# import pandas as pd
import collections
from collections import OrderedDict



def checknormality(data: list[float]):
    stat, p = scipy.stats.shapiro(data)
    print("After calculating p value ", p)
    print('stat=%.3f, p=%.3f' % (stat, p))
    isnormal: bool
    if p > 0.01:
        print('Probably Gaussian')
        isnormal = True
    else:
        print('Probably not Gaussian')
        isnormal = False
    return p, isnormal


def readfeaturesmetadatafile(featuresfilename: str) -> tuple[dict, dict, dict, dict, dict, dict, dict, dict]:
    featurefileheaderlist: list[str] = []
    featurecounter: int
    featurecounter = 0
    featureinfodict: dict = {}
    transcript_or_feature_to_gene_association: dict = {}
    transcript_or_feature_to_confidence_association: dict = {}
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
                    # row[0] is the feature id
                    # header is _id	confidence	count	Pathway	module	Gene Transcript
                    # row[0] is transcript information
                    if row[0] not in featureinfodict.keys():
                        featureinfodict[row[0]] = {}
                    if featurefileheaderlist[counter] not in featureinfodict[row[0]].keys():
                        featureinfodict[row[0]][featurefileheaderlist[counter]] = row[1:][counter]
                    # add the transcript - gene association
                    if featurefileheaderlist[counter] == 'Gene':
                        transcript_or_feature_to_gene_association[row[0]] = row[1:][counter]
                        if row[1:][counter] not in gene_to_transcript_or_feature_association.keys():
                            gene_to_transcript_or_feature_association[row[1:][counter]] = []
                        # Append transcript / feature information
                        gene_to_transcript_or_feature_association[row[1:][counter]].append(row[0])
                    if featurefileheaderlist[counter] == 'confidence':
                        transcript_or_feature_to_confidence_association[row[0]] = row[1:][counter]
                    if featurefileheaderlist[counter] == 'count':
                        transcript_or_feature_to_count_normalized_association[row[0]] = row[1:][counter]
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
        maxcount = max(transcript_or_feature_to_count_normalized_association.values())
        transcript_or_feature_to_count_normalized_association = {k: normalize(v, maxcount) for k, v in
                                                                 transcript_or_feature_to_count_normalized_association.items()}
        import json
        # print(json.dumps(gene_to_transcript_or_feature_association, sort_keys=True, indent=2))
        # print(json.dumps(pathway_to_module_nummodulesperpathway_association, sort_keys=True, indent=2))
        # print(json.dumps(gene_to_pathway_association, sort_keys=True, indent=2))

        for key in pathway_to_gene_association.keys():
            pathway_to_gene_association[key] = list(set(pathway_to_gene_association.get(key)))

        for k, v in pathway_to_module_nummodulesperpathway_association.items():
            print(" Key (Pathway) ", k, " Number of modules ",
                  len(pathway_to_module_nummodulesperpathway_association.get(k)))

        pathway_to_numgenes_normalized: dict = {
            k: normalize(len(v), max(map(len, pathway_to_gene_association.values()))) for k, v in
            pathway_to_gene_association.items()}

        pathway_to_nummodules_normalized: dict = {
            k: normalize(len(v), max(map(len, pathway_to_module_nummodulesperpathway_association.values()))) for k, v in
            pathway_to_module_nummodulesperpathway_association.items()}
        print(json.dumps(pathway_to_numgenes_normalized, sort_keys=True, indent=2))
        print(json.dumps(pathway_to_nummodules_normalized, sort_keys=True, indent=2))
        print(" type of pathay num modules list ", type(pathway_to_nummodules_normalized.values()))

        (pval_pathwaymods, normality_pathwaymods) = checknormality(list(pathway_to_nummodules_normalized.values()))
        (pval_pathwaygenes, normality_pathwaygenes) = checknormality(list(pathway_to_numgenes_normalized.values()))

        # df1 = pd.DataFrame.from_dict(pathway_to_nummodules_normalized)
        # df2 = pd.DataFrame.from_dict(pathway_to_numgenes_normalized)
        # print(df1)
        # print(df2)
        # exit(1)
        if normality_pathwaymods is True and normality_pathwaygenes is True:
            dataforfindingcovariance = list(
                zip(list(pathway_to_nummodules_normalized.values()), list(pathway_to_numgenes_normalized.values())))
            print(dataforfindingcovariance)
            print("After printing zip")
            covariance_between_numberofmodules_in_pathway_to_numberofgenes = numpy.cov(dataforfindingcovariance)
            print(covariance_between_numberofmodules_in_pathway_to_numberofgenes)

            pathwaycounter: int
            pathwaycounter = 0
            pathwayranks_name_to_weighted_gene_correlation_sum: dict[str, float] = {}
            weighted_gene_correlation_vals: float = 0
            for rows in range(len(covariance_between_numberofmodules_in_pathway_to_numberofgenes)):
                pathwayranks_name_to_weighted_gene_correlation_sum[
                    list(pathway_to_numgenes_normalized.keys())[pathwaycounter]] = 0
                weighted_gene_correlation_vals = 0
                for cols in range(len(covariance_between_numberofmodules_in_pathway_to_numberofgenes[rows])):
                    weighted_gene_correlation_vals += \
                        covariance_between_numberofmodules_in_pathway_to_numberofgenes[rows][cols]
                pathwayranks_name_to_weighted_gene_correlation_sum[
                    list(pathway_to_numgenes_normalized.keys())[pathwaycounter]] = weighted_gene_correlation_vals
                pathwaycounter = pathwaycounter + 1

            print(pathwayranks_name_to_weighted_gene_correlation_sum)
            import operator
            # sort in ascending order
            pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes = dict(sorted(pathwayranks_name_to_weighted_gene_correlation_sum.items(), key=lambda item: item[1]))
            print('Dictionary in descending order by value : ',
                  pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes)
            print(type(pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes))

    return featureinfodict, transcript_or_feature_to_gene_association, transcript_or_feature_to_confidence_association, \
           transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
           pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, gene_to_pathway_association, pathway_to_gene_association
