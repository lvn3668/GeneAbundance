# Author Lalitha Viswanathan
# Function to rank genes based on pathways


def rankgenesbasedonpathways(geneid_expressionvalues: dict, gene_to_pathway_association: dict,
                             pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes: dict) -> dict:
    """
    :param geneid_expressionvalues:
    :param gene_to_pathway_association:
    :param pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes:
    :return:
    """
    genesrankedbypathwayexpression: dict = {}
    for gene in geneid_expressionvalues.keys():
        # get rank of the pathways it is expressed in
        if gene not in genesrankedbypathwayexpression.keys():
            genesrankedbypathwayexpression[gene] = 0
        pathwayrank: float = 0
        for pathway in set(gene_to_pathway_association.get(gene)):
            print("Pathway associated with gene *", pathway.strip(), "* ",
              pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes.get(pathway.strip()))
            pathwayrank += pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes.get(pathway.strip())
        genesrankedbypathwayexpression[gene] = pathwayrank
    return genesrankedbypathwayexpression