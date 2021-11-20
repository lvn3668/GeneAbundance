# Author: Lalitha Viswanathan
# Correlate samples to genes

def find_sample_gene_association(transcript_or_feature_to_gene_association: dict,
                                 sample_to_expressed_feature_or_transcript_association: dict):

    sample_to_gene_association: dict = {}
    gene_to_sample_association: dict = {}

    for sampleid in sample_to_expressed_feature_or_transcript_association.keys():
        # get transcripts for a given sample (transcripts and features used interchangeably)
        #print("Sample id ", sampleid.strip())
        #print(" Transcripts ", sample_to_expressed_feature_or_transcript_association.get(sampleid.strip()))
        if sample_to_expressed_feature_or_transcript_association.get(sampleid.strip()) is not None:
            for transcript in sample_to_expressed_feature_or_transcript_association.get(sampleid.strip()):
                # for each transcript, find the associated gene (1: 1 association)
                # sample to feature/transcript -> 1: many
                # feature / transcript : gene -> 1: 1
                # sample to gene: many to many
                #print("Associated genes ", transcript_or_feature_to_gene_association.get(transcript),  " transcript id ", transcript)
                #print(" Gene Id ", transcript_or_feature_to_gene_association.get(transcript), " sample id ", sampleid)
                if sampleid.strip() not in sample_to_gene_association.keys():
                    sample_to_gene_association[sampleid.strip()] = []
                    sample_to_gene_association[sampleid.strip()].append(transcript_or_feature_to_gene_association.get(transcript))
                if transcript_or_feature_to_gene_association.get(transcript) not in gene_to_sample_association.keys():
                    gene_to_sample_association[transcript_or_feature_to_gene_association.get(transcript)] = []
                    gene_to_sample_association[transcript_or_feature_to_gene_association.get(transcript)].append(sampleid.strip())

    import json
    #print(json.dumps(sample_to_gene_association, sort_keys=True, indent=2))
    #print(json.dumps(gene_to_sample_association, sort_keys=True, indent=2))
    return sample_to_gene_association, gene_to_sample_association
