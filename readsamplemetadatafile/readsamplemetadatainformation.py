#Author: Lalitha Viswanathan
# Gene Abundance Detection

# Module to read sample information (metadata information about sample) such as
# barcode sequence, primer sequence, drug-response (yes / no) and days since included in experriment pool

from typing import List, Any, Generator

from readdatafile.readdatafile import normalize


def readsampleinformation(samplesfilename: str, datainfodict2d: dict):
    ###########################################################################
    samplefileheaderlist: list[str] = []
    featurecounter: int
    featurecounter = 0
    missingsamples: list[str] = []
    subjectstreatedwithdrugs: dict = {}
    samplid_to_subject_association: dict = {}
    sampleid_to_dayssinceexperimentstarted: dict = {}
    with open(samplesfilename) as samplefile:
        samplefilerows: Generator[list[str], Any, None] = (line.split('\t') for line in samplefile)
        for row in samplefilerows:
            if featurecounter != 0:
                loopingcounter: int
                for loopingcounter in range(1, len(samplefileheaderlist)):
                    # row[0] is header info
                    if row[0] in datainfodict2d.keys():
                        datainfodict2d[row[0]][samplefileheaderlist[loopingcounter]] = row[1:][loopingcounter]
                        if samplefileheaderlist[loopingcounter] == "Subject":
                            # key into hash is subject1 / subject2 and
                            # treated with drugs yes /no is the 2nd key into the hash
                            # value of this 2d key is sample  id
                            if row[1:][loopingcounter] not in subjectstreatedwithdrugs.keys():
                                subjectstreatedwithdrugs[row[1:][loopingcounter]] = {}
                            # 1D Hash created with subject1 subject2
                            if row[1:][loopingcounter + 1] not in subjectstreatedwithdrugs[
                                row[1:][loopingcounter]].keys():
                                subjectstreatedwithdrugs[row[1:][loopingcounter]][row[1:][loopingcounter + 1]] = []
                            # 2d hash created with yes/no value for treated-with-drugs field
                            subjectstreatedwithdrugs[row[1:][loopingcounter]][row[1:][loopingcounter + 1]].append(
                                row[0])
                            # create association of sample to days since experiment started
                            sampleid_to_dayssinceexperimentstarted[row[0]] = row[1:][loopingcounter+2]
                            # create association of sample id and type of subject (subject 1 / subject 2)
                            if row[0] not in samplid_to_subject_association.keys():
                                samplid_to_subject_association[row[0]] = row[1:][loopingcounter]

                    else:
                        missingsamples.append(row[0])
            elif featurecounter == 0:
                # header row
                for field in row[1:]:
                    # fields called sample, values, barcode sequence, linker primer sequence
                    samplefileheaderlist.append(field)
            featurecounter = featurecounter + 1

    print(len(sampleid_to_dayssinceexperimentstarted.values()))
    maxcount = max(sampleid_to_dayssinceexperimentstarted.values())
    sampleid_to_dayssinceexperimentstarted_normalized = {k: normalize(v, maxcount) for k, v in sampleid_to_dayssinceexperimentstarted.items()}
    return datainfodict2d, samplid_to_subject_association, subjectstreatedwithdrugs, missingsamples,sampleid_to_dayssinceexperimentstarted_normalized
##################################################################################################
