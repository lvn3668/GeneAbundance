# Author Lalitha Viswanathan
# Read fna file
# Gene Abundance detection (Find most abundantly expressed genes)
from readdatafile.readdatafile import normalize


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
    maxcount = max(sampleid_to_length_association.values())
    sampleid_to_length_association_normalizedovermaxlength = {k: normalize(v, maxcount) for k, v in sampleid_to_length_association.items()}
    return sequencedict, sampleid_to_length_association_normalizedovermaxlength
###########################################################################