def findrelativeabundanceofgenes(subjectstreatedwithdrugs: dict, expressedtranscriptspersample: dict,
                                 transcript_gene_assoc: dict, relativeabundancepersample_sorted: dict,
                                 genesrelativeabundance: dict):
    geneswithmostrelativeabundance: dict = {}
    for subject in subjectstreatedwithdrugs.keys():
        if subject not in geneswithmostrelativeabundance.keys():
            geneswithmostrelativeabundance[subject] = {}
        treatedwithdrugsyesno: str
        for treatedwithdrugsyesno in subjectstreatedwithdrugs.get(subject).keys():
            if treatedwithdrugsyesno not in geneswithmostrelativeabundance[subject].keys():
                geneswithmostrelativeabundance[subject][treatedwithdrugsyesno] = {}
            for sampleid in subjectstreatedwithdrugs.get(subject).get(treatedwithdrugsyesno):
                if sampleid not in geneswithmostrelativeabundance[subject][treatedwithdrugsyesno].keys():
                    geneswithmostrelativeabundance[subject][treatedwithdrugsyesno][sampleid] = {}
                if sampleid in expressedtranscriptspersample.keys():
                    for transcript in expressedtranscriptspersample.get(sampleid):
                        if transcript not in geneswithmostrelativeabundance[subject][treatedwithdrugsyesno][sampleid].keys():
                            geneswithmostrelativeabundance[subject][treatedwithdrugsyesno][sampleid][transcript] = {}
                            if transcript_gene_assoc.get(transcript) not in \
                                    geneswithmostrelativeabundance[subject][treatedwithdrugsyesno][sampleid].keys():
                                geneswithmostrelativeabundance[subject][treatedwithdrugsyesno][sampleid][transcript][
                                    transcript_gene_assoc.get(transcript)] = {}
                                geneswithmostrelativeabundance[subject][treatedwithdrugsyesno][sampleid][transcript][
                                    transcript_gene_assoc.get(transcript)] = relativeabundancepersample_sorted.get(sampleid)
                                if transcript_gene_assoc.get(transcript) not in genesrelativeabundance.keys():
                                    genesrelativeabundance[transcript_gene_assoc.get(transcript)] = []
                                    # Add mean relative abundance to generelativeabundance
                                    # gene expressed via multiple transcripts in multiple samples
                                    # hence list
                                genesrelativeabundance[transcript_gene_assoc.get(transcript)].append(
                                    relativeabundancepersample_sorted.get(sampleid))
    return genesrelativeabundance