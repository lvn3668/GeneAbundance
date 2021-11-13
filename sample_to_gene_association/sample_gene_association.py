# Author: Lalitha Viswanathan
# Correlate samples to genes

def remove_duplicates(dictvaluewithduplicates):
    return list(set(dictvaluewithduplicates))

def find_sample_gene_association(transcript_or_feature_to_gene_association: dict, sample_to_expressed_feature_or_transcript_association: dict):
    # key is transcript or feature
    sample_to_gene_association: dict = {}
    gene_to_sample_association: dict = {}
    for sampleid in sample_to_expressed_feature_or_transcript_association.keys():
        # get transcripts for a given sample (transcripts and features used interchangeably)
        for key in sample_to_expressed_feature_or_transcript_association.get(sampleid):
            # for each transcript, find the associated gene (1: 1 association)
            # sample to feature/transcript -> 1: many
            # feature / transcript : gene -> 1: 1
            # sample to gene: many to many
            if type(transcript_or_feature_to_gene_association.get(key)) == "<class 'list'>":
                for geneid in list(transcript_or_feature_to_gene_association.get(key)):
                    if sampleid not in sample_to_gene_association.keys():
                        sample_to_gene_association[sampleid] = []
                    sample_to_gene_association[sampleid].append(geneid)
                    if geneid not in gene_to_sample_association.keys():
                        gene_to_sample_association[geneid] = []
                    gene_to_sample_association[geneid].append(sampleid)
            else:
                if sampleid not in sample_to_gene_association.keys():
                        sample_to_gene_association[sampleid] = transcript_or_feature_to_gene_association.get(key)
                if transcript_or_feature_to_gene_association.get(key) not in gene_to_sample_association.keys():
                        gene_to_sample_association[transcript_or_feature_to_gene_association.get(key)] = sampleid

    import json
    print(json.dumps(sample_to_gene_association, sort_keys=True, indent=2))
    print(json.dumps(gene_to_sample_association, sort_keys=True, indent=2))

    return sample_to_gene_association, gene_to_sample_association