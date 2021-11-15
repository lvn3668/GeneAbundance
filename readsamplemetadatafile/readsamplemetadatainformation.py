#Author: Lalitha Viswanathan
# Gene Abundance Detection

# Module to read sample information (metadata information about sample) such as
# barcode sequence, primer sequence, drug-response (yes / no) and days since included in experriment pool
import json
from typing import Any, Generator

from datanormalization.normalizeovermaxvalue import normalizeovermaximum

# Read sample metadata
def readsamplemetadatainformation(samplesfilename: str,
                                  sample_expressedfeature_or_transcript_expressionvalues_matrix: dict) -> tuple[
    dict, dict, dict, list[str], dict]:
    """

    :rtype: object
    :param samplesfilename:
    :param sample_expressedfeature_or_transcript_expressionvalues_matrix:
    :return:
    """
    ###########################################################################
    samplefileheaderlist: list[str] = []
    samplerowcounter: int
    samplerowcounter = 0
    samples_with_missing_expression_values: list[str] = []
    subjects_treatedwithdrugs_to_sample_association: dict = {}
    samplid_to_subject_association: dict = {}
    sampleid_to_dayssinceexperimentstarted_association: dict = {}
    with open(samplesfilename) as samplefile:
        samplemetadatarows: Generator[list[str], Any, None] = (line.split('\t') for line in samplefile)
        for persamplemetadatainformation in samplemetadatarows:
            if samplerowcounter != 0:
                counter: int
                for counter in range(1, len(samplefileheaderlist)):
                    # persamplemetadatainformation[0] is header info
                    if persamplemetadatainformation[0] in sample_expressedfeature_or_transcript_expressionvalues_matrix.keys():
                        # first index is the sample id (which is the key when the expression values file was read in)
                        # 2nd index is the header information such as barcode sequence, linker sequence
                        sample_expressedfeature_or_transcript_expressionvalues_matrix[persamplemetadatainformation[0]][samplefileheaderlist[counter]] = persamplemetadatainformation[1:][counter]
                        if samplefileheaderlist[counter] == "Subject":
                            # key into hash is subject1 / subject2 and
                            # treated with drugs yes /no is the 2nd key into the hash
                            # value of this 2d key is sample  id
                            if persamplemetadatainformation[1:][counter] not in subjects_treatedwithdrugs_to_sample_association.keys():
                                subjects_treatedwithdrugs_to_sample_association[persamplemetadatainformation[1:][counter]] = {}
                            # 1D Hash created with subject1 subject2
                            if persamplemetadatainformation[1:][counter + 1] not in subjects_treatedwithdrugs_to_sample_association[
                                persamplemetadatainformation[1:][counter]].keys():
                                subjects_treatedwithdrugs_to_sample_association[persamplemetadatainformation[1:][counter]][persamplemetadatainformation[1:][counter + 1]] = []
                            # 2d hash created with yes/no value for treated-with-drugs field
                            subjects_treatedwithdrugs_to_sample_association[persamplemetadatainformation[1:][counter]][persamplemetadatainformation[1:][counter + 1]].append(
                                persamplemetadatainformation[0])
                            # create association of sample to days since experiment started
                            sampleid_to_dayssinceexperimentstarted_association[persamplemetadatainformation[0]] = persamplemetadatainformation[1:][counter+2]
                            # create association of sample id and type of subject (subject 1 / subject 2)
                            if persamplemetadatainformation[0] not in samplid_to_subject_association.keys():
                                samplid_to_subject_association[persamplemetadatainformation[0]] = persamplemetadatainformation[1:][counter]
                    else:
                        samples_with_missing_expression_values.append(persamplemetadatainformation[0])
            elif samplerowcounter == 0:
                # header persamplemetadatainformation
                for field in persamplemetadatainformation[1:]:
                    # fields called sample, values, barcode sequence, linker primer sequence
                    samplefileheaderlist.append(field)
            samplerowcounter = samplerowcounter + 1

    print(len(sampleid_to_dayssinceexperimentstarted_association.values()))
    sampleid_to_dayssinceexperimentstarted_normalized = normalizeovermaximum(sampleid_to_dayssinceexperimentstarted_association)

    return sample_expressedfeature_or_transcript_expressionvalues_matrix, samplid_to_subject_association, \
           subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
           sampleid_to_dayssinceexperimentstarted_normalized
##################################################################################################
