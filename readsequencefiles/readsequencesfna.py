# Author Lalitha Viswanathan
# Read fna file
# Gene Abundance detection (Find most abundantly expressed genes)
from datanormalization.normalizeovermaxvalue import normalize_over_maximum_value
from readdatafile.readdatafile import normalize

# Author: Lalitha Viswanathan
# Read sequences from fna File
# and normalize sequence length over max length
def readfnafile(sequencesfilename: str):
    sequencedict: dict = {}
    header: str
    sampleid_to_length_association: dict = {}
    with open(sequencesfilename) as file:
        while line := file.readline().rstrip():
            if line.startswith(">"):
                header = line
                sampleid = header.split(r'\s')[0].replace(">", "").split("_")[0]
                if sampleid not in sampleid_to_length_association.keys():
                    sampleid_to_length_association[sampleid] = ""
            else:
                sequencedict[header] = line
                sampleid_to_length_association[sampleid] = len(line)
    sampleid_to_length_association_normalizedovermaxlength = normalize_over_maximum_value(sampleid_to_length_association)
    return sequencedict, sampleid_to_length_association_normalizedovermaxlength
###########################################################################