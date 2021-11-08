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
from findrelativeabundanceofgenes.findrelabundancegenes import findrelativeabundanceofgenes
from readdatafile.readdatafile import readdatafile
from readfeaturefile.readfeaturesfiles import readfeaturesfile
from readsamplemetadatafile.readsamplemetadatainformation import readsampleinformation
from readsequencefiles.readsequencesfna import readfnafile
from reversecomplementDNA.reversecomplementDNA import reversecomplementDNA, validationofreversecomplement

# Explore the `samples` metadata. Calculate how many sequenced samples are available
# for each subject (\"Subject\" column) across treatment usage (\"Treated_with_drug\" column).
# Hint: You can use a \"groupby\" function or a \"table\" function.**"

# Transform the `data` table to relative abundances - this is done by dividing the counts in each sample and
# each feature by the sum of the counts for all features in each sample. Name this new dataframe `data_ra`.**"
# Determine the top 10 genes with the highest mean relative abundance found across both subjects.
# The Gene label can be found in the `features` dataframe.**"

def findgenesofsignificance(datafilename: str, featuresfilename: str, samplesfilename: str, sequencesfilename: str):

    sequencedict, sampleid_to_length_association = readfnafile(sequencesfilename)
    geneexprvalsperfeaturespersample, relativeabundancepersample, expressedtranscriptspersample = readdatafile(datafilename)
    #         datainfodict2d, samplid_to_subject_association, subjectstreatedwithdrugs, missingsamples, sampleid_to_dayssinceexperimentstarted_normalized
    geneexprvalsperfeaturespersample, samplessubject, subjectstreatedwithdrugs, missingsamples, sampleid_to_dayssinceexperimentstarted_normalized = readsampleinformation(samplesfilename, geneexprvalsperfeaturespersample)
    featureinfodict2d, transcriptgene,transcriptconfidence, transcriptcount  = readfeaturesfile(featuresfilename)
    relativeabundancepersample_sorted = dict(sorted(relativeabundancepersample.items(),
                                                    key=lambda item: item[1],
                                                    reverse=True))
    genesrelativeabundance: dict = {}
    print("######### relative abundance per sample keys ")
    print(relativeabundancepersample_sorted.keys())
    print("######### transcript gene ")
    print(transcriptgene.keys())
    print("##### subject treated with drugs ")
    print(subjectstreatedwithdrugs.keys())
    genesrelativeabundance = findrelativeabundanceofgenes(subjectstreatedwithdrugs, expressedtranscriptspersample,
                                 transcriptgene, relativeabundancepersample_sorted,
                                 genesrelativeabundance)

    intersection = dict(relativeabundancepersample_sorted.items() & transcriptgene.items())
    print(intersection)

    for gene in genesrelativeabundance.keys():
        # Multiply by confidence to normalize
        # Use counts to normalize
        listofabundancevals: list[str] = []
        for relabundance in genesrelativeabundance.get(gene):
            if gene in transcriptconfidence.keys():
                relabundance = relabundance * float(transcriptconfidence.get(gene))
            if gene in transcriptcount.keys():
                relabundance = relabundance * float(transcriptcount.get(gene))
        # abundance values normalized by count and confidence level
        # Use days since experiment started for treated-with-drug NO
        # If it is NOT treated with drugs and days since experiment started is greater than 0
            # Inverse proportion
            # Longer the ample in the experiment, and not treated with drug, over expressed it is,
            # Normalize days since experiment started (0/ (max), 10 /(max) , etc.)
            # Abundance is unaltered for treated-with-drug yes cases (????)
            # For untreated cases, account for over expression / underexpression by DIVIDING the abundance measure
            # by normalized days-since-experiment-started
            # This will underexpress those transcripts belonging to UNTREATED samples, that have been the longest in the experiment pool

            if gene in sampleid_to_dayssinceexperimentstarted_normalized.keys() and sampleid_to_dayssinceexperimentstarted_normalized.get(gene) != 0:
                relabundance = relabundance / sampleid_to_dayssinceexperimentstarted_normalized.get(gene)
                listofabundancevals.append(relabundance)
        # Sample length from fna file
        # check for correlation between sample length and number of expressed transcripts / features
        # assuming positive correlation


        # normalize over sample length (gene abundance can come from multiple transcripts from different samples)

        # build pathway module gene association
        # A gene can be implicated in multiple pathways and a pathway should have mutiple modules
        # A gene from a well expressed pathway -> HOW to normalize
        genesrelativeabundance[gene] = listofabundancevals
        print("Gene ", gene, " Relative abundances ", genesrelativeabundance.get(gene))

    # Determine the top 10 genes with the highest mean relative abundance found across both subjects.
    # The Gene label can be found in the `features` dataframe.**"

    ##################################################################################################
    return sequencedict, featureinfodict2d, geneexprvalsperfeaturespersample, missingsamples, subjectstreatedwithdrugs, relativeabundancepersample_sorted


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
                print(line)
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
        sequencedict, featureinfodict2d, datainfodict2d, missingsamples, subjectstreatedwithdrugs, relativeabundancepersample = findgenesofsignificance(
            datafilename, featuresfilename, samplesfilename, sequencesfilename)
        print("After reading in samples, data features")
        for key in datainfodict2d.keys():
            print(key)
            print(datainfodict2d.get(key))
        print("Relative abundance printing ")
        for key in relativeabundancepersample:
            print("***", key, "\t", relativeabundancepersample.get(key))
    finally:
        print(list(checkforvalidparentheses("{([((()])))[][[()]]}")))
        print("After printing the first 15 lines")
