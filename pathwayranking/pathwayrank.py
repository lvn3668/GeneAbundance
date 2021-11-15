# Author: Lalitha Viswanathan
# Pathway Ranking by weighted correlation of number of modules in each pathway and number of genes in each pathway
import json
from typing import Any
import numpy

from checkforgaussiandistribution.gaussiandistributionusingshapirowilkis import checkforgaussiandistribution
from readdatafile.readdatafile import normalize

# Return ranking of pathways in features.txt file in descending order
def pathwayrank(pathway_to_gene_association: dict, pathway_to_module_nummodulesperpathway_association: dict
                ) -> dict:
    pathway_to_numgenes_normalized: dict
    pathway_to_nummodules_normalized: dict
    for key in pathway_to_gene_association.keys():
        pathway_to_gene_association[key] = list(set(pathway_to_gene_association.get(key)))
    pathway_to_numgenes_normalized: dict = {
        k: normalize(len(v), max(map(len, pathway_to_gene_association.values()))) for k, v in
        pathway_to_gene_association.items()}
    ###################################################################################################
    # Normalize number of modules per pathway
    for k, v in pathway_to_module_nummodulesperpathway_association.items():
        print(" Key (Pathway) ", k, " Number of modules ",
              len(pathway_to_module_nummodulesperpathway_association.get(k)))
    pathway_to_nummodules_normalized: dict = {
        k: normalize(len(v), max(map(len, pathway_to_module_nummodulesperpathway_association.values()))) for k, v in
        pathway_to_module_nummodulesperpathway_association.items()}

    ###################################################################################################
    print(json.dumps(pathway_to_numgenes_normalized, sort_keys=True, indent=2))
    print(json.dumps(pathway_to_nummodules_normalized, sort_keys=True, indent=2))
    print(" type of pathay num modules list ", type(pathway_to_nummodules_normalized.values()))

    ###################################################################################################
    # Check for normality of data
    (pval_pathwaymods, normality_pathwaymods) = checkforgaussiandistribution(
        list(pathway_to_nummodules_normalized.values()))
    (pval_pathwaygenes, normality_pathwaygenes) = checkforgaussiandistribution(
        list(pathway_to_numgenes_normalized.values()))
    # Find correlation between number of modules per pathway and number of genes per pathway
    if normality_pathwaymods is True and normality_pathwaygenes is True:
        dataforfindingcovariance: list[tuple[Any, Any]] = list(
            zip(list(pathway_to_nummodules_normalized.values()), list(pathway_to_numgenes_normalized.values())))
        print(dataforfindingcovariance)
        covariance_between_numberofmodules_in_pathway_to_numberofgenes = numpy.cov(dataforfindingcovariance)
        print(covariance_between_numberofmodules_in_pathway_to_numberofgenes)
        # loop through the covariance matrix and add the weights
        # This weighted sum is presently used to rank pathways
        pathwaycounter: int
        pathwaycounter = 0
        pathwayranks_name_to_weighted_gene_correlation_sum: dict[str, float] = {}
        weighted_gene_correlation_vals: float
        for rows in range(len(covariance_between_numberofmodules_in_pathway_to_numberofgenes)):
            pathwayranks_name_to_weighted_gene_correlation_sum[
                list(pathway_to_numgenes_normalized.keys())[pathwaycounter]] = 0
            weighted_gene_correlation_vals = 0
            for cols in range(len(covariance_between_numberofmodules_in_pathway_to_numberofgenes[rows])):
                weighted_gene_correlation_vals += \
                    covariance_between_numberofmodules_in_pathway_to_numberofgenes[rows][cols]
            # Dict of pathway Id and weighted gene correlation vals
            pathwayranks_name_to_weighted_gene_correlation_sum[
                list(pathway_to_numgenes_normalized.keys())[pathwaycounter]] = weighted_gene_correlation_vals
            pathwaycounter = pathwaycounter + 1

        ###################################################################################################

        # Rank pathways by weighted sum
        pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes = dict(sorted(pathwayranks_name_to_weighted_gene_correlation_sum.items(), key=lambda item: item[1]))
        print(pathwayranks_name_to_weighted_gene_correlation_sum)
        print('Dictionary in descending order by value : ',
              pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes)
        print(type(pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes))
        ###################################################################################################
    return pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes
