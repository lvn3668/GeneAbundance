# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

# "Write a function to evaluate whether parentheses are matching across a string.
# In other words each open bracket has a matching closing bracket and are correctly ordered.
# The input should be a string and the return should be a boolean - `True` if the brackets are correct `False` if they are not.
# The only characters used will be `'{}[]()'`"


# `samples` and `features` tables contain metadata about the rows and columns of the `data` table.
# Follow the instructions in your preferred programming language to answer the following questions."
# Create a variable named `data` by reading in a dataframe from the file `data.txt`
# which contains the sample (col) by feature counts (rows). The first column named `_id` can be used as index.**"
# Create a variable named `samples` by reading in a dataframe from the file `samples.txt` with the sample metadata.
# Add a checkpoint to verify that the `#SampleID` column has unique values.**"
# Create a variable named `features` by reading in a dataframe from the file `features.txt` with the features metadata.
# The first column named `_id` can be used as index.**"

import argparse
import re
# "Imagine a grocery store with multiple aisles and multiple people waiting in each line.
# You quickly look at the items in each customers basket and know exactly how long it will take them to check out.
# If each customer always chooses the shortest line, calculate how long it will take for all customers to go through checkout.
# For example if you have 3 customers who will take 5, 10, and 3 minutes each to go through checkout and there is only one lane available it will take 18 minutes,
# if there are two aisles available it will take ten minutes"
# Write a method to reverse translate the string “ACCTGGCCGTACCT”.
# Results should be “AGGTACGGCCAGGT”.
# Add a validation step to show that the results is what expected.\n",
from checkfornumberofperfectsquares.countperfectsquares import findperfectsquares
from findexpressionvalues_sample_gene_association.findexpressionvalues_sample_gene_association import \
    findexpressionvalues_gene_sample_association
from readdatafile.readdatafile import readsampletofeatureexpressionvalues, normalize
from readfeaturefile.readfeaturesfiles import readfeaturesmetadatafile
from readsamplemetadatafile.readsamplemetadatainformation import readsamplemetadatainformation
from readsequencefiles.readsequencesfna import readfnafile
from reversecomplementDNA.reversecomplementDNA import reversecomplementDNA, validationofreversecomplement
from sample_to_gene_association.sample_gene_association import *


# Explore the `samples` metadata. Calculate how many sequenced samples are available
# for each subject (\"Subject\" column) across treatment usage (\"Treated_with_drug\" column).
# Hint: You can use a \"groupby\" function or a \"table\" function.**"

# Transform the `data` table to relative abundances - this is done by dividing the counts in each sample and
# each feature by the sum of the counts for all features in each sample. Name this new dataframe `data_ra`.**"
# Determine the top 10 genes with the highest mean relative abundance found across both subjects.
# The Gene label can be found in the `features` dataframe.**"

def findgenesofsignificance(datafilename: str, featuresfilename: str, samplesfilename: str, sequencesfilename: str):
    print("Inside find genes of significance ", sequencesfilename)
    sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength = readfnafile(
        sequencesfilename)

    sample_expressedfeature_or_transcript_expressionvalues_matrix, meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
    expressedtranscriptspersample_sampletotranscriptassoc = readsampletofeatureexpressionvalues(
        datafilename)

    sample_expressedfeature_or_transcript_expressionvalues_matrix, samplid_to_subject_association, \
    subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
    sampleid_to_dayssinceexperimentstarted_normalized = readsamplemetadatainformation(
        samplesfilename, sample_expressedfeature_or_transcript_expressionvalues_matrix)

    features_metadata, transcript_or_feature_to_gene_association, transcript_or_feature_to_confidence_association, \
    transcript_or_feature_to_count_association, gene_to_transcript_or_feature_association, \
    pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes,\
        gene_to_pathway_association, pathway_to_gene_association = readfeaturesmetadatafile(
        featuresfilename)

    print("Before finding sample to gene association")
    sample_to_gene_association, gene_to_sample_association = find_sample_gene_association(
        transcript_or_feature_to_gene_association, expressedtranscriptspersample_sampletotranscriptassoc)

    print("######### transcript gene ")
    print(transcript_or_feature_to_gene_association.keys())
    print("##### subject treated with drugs ")
    print(subjects_treatedwithdrugs_to_sample_association.keys())
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene: dict
    geneid_expressionvalues: dict
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues \
        = findexpressionvalues_gene_sample_association(
        subjects_treatedwithdrugs_to_sample_association,
        expressedtranscriptspersample_sampletotranscriptassoc,
        transcript_or_feature_to_gene_association,
        meanexpression_per_sample_for_all_features_or_transcripts_sorted)

    intersection = dict(
        meanexpression_per_sample_for_all_features_or_transcripts_sorted.items() & transcript_or_feature_to_gene_association.items())
    print(intersection)

    print('Dictionary in descending order by value : ',
          pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes)

    genesrankedbypathwayexpression: dict = {}
    for gene in geneid_expressionvalues.keys():
        # get rank of the pathways it is expressed in
        if gene not in genesrankedbypathwayexpression.keys():
            genesrankedbypathwayexpression[gene] = 0

        pathwayrank: float = 0
        print("Pathway associated with gene ", gene, " ", gene_to_pathway_association.get(gene))
        for pathway in set(gene_to_pathway_association.get(gene)):
            print("Pathway associated with gene *", pathway.strip(),"* ", pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes.get(pathway.strip()))
            pathwayrank += pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes.get(pathway.strip())
        print("Gene ", gene, " Pathway rank ", pathwayrank)
        genesrankedbypathwayexpression[gene] = pathwayrank

    import json
    print(json.dumps(genesrankedbypathwayexpression, sort_keys=True, indent=2))
    print("After printing gene list ranked by pathway expression NOT normalized ( num genes num modules correlated)")

    genesmissingpathwayassociations: list[str] = []
    for gene in geneid_expressionvalues.keys():
        # Multiply by confidence to normalize
        # Use counts to normalize
        listofabundancevals: list[str] = []
        for expressionvaluesnormalizedoverexpressedfeatures_or_transcripts in geneid_expressionvalues.get(gene):
            if genesrankedbypathwayexpression.get(gene) != 0:
                expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts*genesrankedbypathwayexpression.get(gene)
            else:
                genesmissingpathwayassociations.append(gene)
            # key is gene
            # value is normalized expression values
            for transcript in gene_to_transcript_or_feature_association.get(gene):
                # sum all the transcript confidence measures
                sum_transcript_confidence_measure: float = 0
                if transcript in transcript_or_feature_to_confidence_association.keys():
                    sum_transcript_confidence_measure = sum_transcript_confidence_measure + float(transcript_or_feature_to_confidence_association.get(transcript))

                # then normalize by count
                sum_transcript_count_measure: float = 0
                if transcript in transcript_or_feature_to_count_association.keys():
                    sum_transcript_count_measure = sum_transcript_count_measure + float(transcript_or_feature_to_count_association.get(transcript))

            # normalize by confidence and count (positive correlation)
            expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts * sum_transcript_confidence_measure * sum_transcript_count_measure

            # abundance values normalized by count and confidence level
            # Use days since experiment started for treated-with-drug NO
            # If it is NOT treated with drugs and days since experiment started is greater than 0
            # Inverse proportion
            # Longer the sample in the experiment, and not treated with drug, over expressed it is,
            # Normalize days since experiment started (0/ (max), 10 /(max) , etc.)
            # Abundance is unaltered for treated-with-drug yes cases (????)
            # For untreated cases, account for over expression / under expression by DIVIDING the abundance measure
            # by normalized days-since-experiment-started
            # This will underexpress those transcripts belonging to UNTREATED samples,
            # that have been the longest in the experiment pool
            for sample in gene_to_sample_association.get(gene):
                sumdayssinceexperimentstarted: float = 0
                if sample in sampleid_to_dayssinceexperimentstarted_normalized.keys() and \
                        sampleid_to_dayssinceexperimentstarted_normalized.get(sample) != 0:
                    sumdayssinceexperimentstarted = sumdayssinceexperimentstarted + sampleid_to_dayssinceexperimentstarted_normalized.get(sample)

            if sumdayssinceexperimentstarted !=0:
                expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts / sumdayssinceexperimentstarted

            # gene is part of multiple samples
            # for each sample, that the gene is expressed in
            # get sequence length (normalized over longest sequence length)
            # aggregate normalized sequence lengths
            # divide by the sum
            # negative correlation (longer the sequence, longer the PCR fragments and higher expression values)
            # correction is by dividing the sum of lengths
            for sample in gene_to_sample_association.get(gene):
                    sumnormalizedsequencelengths: float = 0
                    if sample in sampleid_to_length_association_normalizedoverrmaxlength.keys() and \
                            sampleid_to_length_association_normalizedoverrmaxlength.get(sample) != 0:
                        sumnormalizedsequencelengths = sumnormalizedsequencelengths + sampleid_to_length_association_normalizedoverrmaxlength.get(
                            sample)

            if sumnormalizedsequencelengths !=0:
                expressionvaluesnormalizedoverexpressedfeatures_or_transcripts = expressionvaluesnormalizedoverexpressedfeatures_or_transcripts / sumnormalizedsequencelengths
            listofabundancevals.append(expressionvaluesnormalizedoverexpressedfeatures_or_transcripts)



        # TBD:
        # 1. Use numpy and panda to find correlations between
        # sequence length and number of expressed transcripts
        # number of expressed trancripts / features and days in pool
        # ???? anything between sequence length and days in pool (positive correln?)
        # transcript length and confidence measure (positive correln expected)
        # transcript count????
        # correct normaliztion algo based on the correln coefficients
        # Sample length from fna file
        # check for correlation between sample length and number of expressed transcripts / features
        # assuming positive correlation

        # normalize over sample length (gene abundance can come from multiple transcripts from different samples)

        ## Pathway ranking (done)
        # build pathway module gene association
        # A gene can be implicated in multiple pathways and a pathway should have mutiple modules
        # A gene from a well expressed pathway -> HOW to normalize

        # to check: functional association and gene abundance expression
        geneid_expressionvalues[gene] = listofabundancevals
        print("Gene ", gene, " Relative abundances ", geneid_expressionvalues.get(gene))

    # Determine the top 10 genes with the highest mean relative abundance found across both subjects.
    # The Gene label can be found in the `features` dataframe.**"
    ##################################################################################################
    return sequence_header_to_fnasequence_association, features_metadata, sample_expressedfeature_or_transcript_expressionvalues_matrix, samples_with_missing_expression_values, subjects_treatedwithdrugs_to_sample_association


def ParseNestedParen(string, level):
    """
    Generate strings contained in nested (), indexing i = level
    """
    if len(re.findall("\(", string)) == len(re.findall("\)", string)):
        LeftRightIndex = [x for x in zip(
            [Left.start() + 1 for Left in re.finditer('\(', string)],
            reversed([Right.start() for Right in re.finditer('\)', string)]))]

    elif len(re.findall("\(", string)) > len(re.findall("\)", string)):
        return ParseNestedParen(string + ')', level)

    elif len(re.findall("\(", string)) < len(re.findall("\)", string)):
        return ParseNestedParen('(' + string, level)

    else:
        return 'fail'

    return [string[LeftRightIndex[level][0]:LeftRightIndex[level][1]]]


def checkforvalidparentheses(inputstr: str):
    print(inputstr)
    stack = []
    for i, c in enumerate(inputstr):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            start = stack.pop()
            yield len(stack), inputstr[start + 1: i]


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    # open sequences.fna and print the first 15 lines
    # See PyCharm help at https://www.jetbrains.com/help/pycharm/
    counter = 1

    parser = argparse.ArgumentParser(description='Process some integers.')
    # parser.add_argument('revcomplementstring', metavar='N', type=str, nargs='+',
    #                    help='string to be reverse complement')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                    const=sum, default=max,
    #                    help='sum the integers (default: find the max)')

    args = parser.parse_args()
    # print(args.accumulate(args.integers))
    try:
        with open("data/sequences.fna") as file:
            while line := file.readline().rstrip():
                # print(line)
                counter = counter + 1
                if counter == 15:
                    break

        print("Before calling reverse complement ")
        DNAString = "ATGC"
        revcomplement = reversecomplementDNA(DNAString)
        print(validationofreversecomplement(DNAString, revcomplement))
        print("After printing reverse complement validation")
        print("Before finding perfect squares")
        findperfectsquares(4, 9)
        findperfectsquares(100, 1748937)
        # findperfectsquares(1341, 74189027341)
        datafilename = "data/data.txt"
        featuresfilename = "data/features.txt"
        samplesfilename = "data/samples.txt"
        sequencesfilename = "data/sequences.fna"
        print("Before finding genes of significance")
        sequencedict, featureinfodict2d, datainfodict2d, missingsamples, subjectstreatedwithdrugs = findgenesofsignificance(
            datafilename, featuresfilename, samplesfilename, sequencesfilename)
        print("After reading in samples, data features")
        #for key in datainfodict2d.keys():
        #    print(key)
        #    print(datainfodict2d.get(key))
        print("Relative abundance printing ")
        # for key in relativeabundancepersample:
        #    print("***", key, "\t", relativeabundancepersample.get(key))
    finally:
        print(list(checkforvalidparentheses("{([((()])))[][[()]]}")))
        print("After printing the first 15 lines")
