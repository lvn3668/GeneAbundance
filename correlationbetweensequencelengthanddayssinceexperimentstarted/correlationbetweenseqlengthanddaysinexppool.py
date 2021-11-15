
def findcorrelationbetweentrancriptlengthanddaysinpool(gene_to_sample_association, sampleid_to_dayssinceexperimentstarted_normalized,
                                                       sampleid_to_length_association_normalizedoverrmaxlength):
    intersectionofsamplelengthanddaysinpool = dict(
        sampleid_to_length_association_normalizedoverrmaxlength.items() & sampleid_to_dayssinceexperimentstarted_normalized.items())
    print(intersectionofsamplelengthanddaysinpool)
    #for gene in gene_to_sample_association.keys():
    #    intersectionofsamplelengthanddaysinpool = dict(sampleid_to_length_association_normalizedoverrmaxlength.items() & sampleid_to_dayssinceexperimentstarted_normalized.items())


