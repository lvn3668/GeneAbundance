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
from reversecomplementDNA.reversecomplementDNA import reversecomplementDNA, validationofreversecomplement


# Explore the `samples` metadata. Calculate how many sequenced samples are available
# for each subject (\"Subject\" column) across treatment usage (\"Treated_with_drug\" column).
# Hint: You can use a \"groupby\" function or a \"table\" function.**"


# Transform the `data` table to relative abundances - this is done by dividing the counts in each sample and
# each feature by the sum of the counts for all features in each sample. Name this new dataframe `data_ra`.**"
# Determine the top 10 genes with the highest mean relative abundance found across both subjects.
# The Gene label can be found in the `features` dataframe.**"

def findgenesofsignificance(datafilename: str, featuresfilename: str, samplesfilename: str, sequencesfilename: str):
    ###########################################################################
    sequencedict: dict = {}
    header: str
    with open(sequencesfilename) as file:
        while line := file.readline().rstrip():
            if line.startswith(">"):
                header = line
            else:
                sequencedict[header] = line
    print("After reading sequences fna file")
    ###########################################################################
    datainfodict2d: dict = {}
    samplelist: list[str] = []
    samplecounter: int
    samplecounter = 0
    with open(datafilename) as datafile:
        rows = (line.split('\t') for line in datafile)
        for row in rows:
            #print(row)
            if samplecounter != 0:
                loopingcounter: int
                for loopingcounter in range(1, len(samplelist)):
                    # first index is the sample id
                    datainfodict2d[samplelist[loopingcounter]][row[0]] = row[1:][loopingcounter]
            elif samplecounter == 0:
                # header row
                # parse to get sample information
                for sample in row[1:]:
                    samplelist.append(sample)
                    datainfodict2d[sample] = {}
            samplecounter = samplecounter + 1
        print(" Number of data file entries read ", samplecounter)
        ###########################################################################
        samplefileheaderlist: list[str] = []
        samplecounter: int
        samplecounter = 0
        missingsamples: list[str] = []
        subjectstreatedwithdrugs: dict = {}
        with open(samplesfilename) as samplefile:
            rows = (line.split('\t') for line in samplefile)
            for row in rows:
                if samplecounter != 0:
                    loopingcounter: int
                    for loopingcounter in range(1, len(samplefileheaderlist)):
                        # row[0] is sample information  / LS108 etc.
                        if row[0] in datainfodict2d.keys():
                            datainfodict2d[row[0]][samplefileheaderlist[loopingcounter]] = row[1:][loopingcounter]
                            if samplefileheaderlist[loopingcounter] == "Subject":
                                # Contains the subject information
                                # Add the sample information
                                # Add dict of subject -> (sample id, treated with drug, days since experiment started)
                                # "Treated_with_drug": row[1:][loopingcounter+1]
                                # key is subject and treatment with drug (yes or no)
                                # value is sample ids
                                if row[1:][loopingcounter] not in subjectstreatedwithdrugs.keys():
                                    subjectstreatedwithdrugs[row[1:][loopingcounter]] = {}
                                if row[1:][loopingcounter+1] not in subjectstreatedwithdrugs[row[1:][loopingcounter]].keys():
                                    subjectstreatedwithdrugs[row[1:][loopingcounter]][row[1:][loopingcounter+1]] = []

                                subjectstreatedwithdrugs[row[1:][loopingcounter]][row[1:][loopingcounter+1]].append(row[0])
                        else:
                            missingsamples.append(row[0])
                elif samplecounter == 0:
                    # header row
                    for field in row[1:]:
                        # fields called sample, values, barcode sequence, linker primer sequence
                        samplefileheaderlist.append(field)

                samplecounter = samplecounter + 1
    print("After reading in sample info records ", samplecounter)
    ##################################################################################################
    featurefileheaderlist: list[str] = []
    featurecounter: int
    featurecounter = 0
    featureinfodict2d: dict = {}
    # 4	module_6	gene_480	transcript_480.6
    # 4483174	1	3	pathway_5	module_67	gene_471	transcript_471.4
    pattern = re.compile("^(\d+\s+\d+\s+\d+\s+\W+\s+\W+\s+\W+\s+\W+)$")
    with open(featuresfilename) as featurefile:
        # leaves out the line 4 module_6 gene_480 transcript_480.6
        # where confidence and count are unavailable
        rows = (line.split('\t') for line in featurefile if pattern.match(line))
        for row in rows:
            if featurecounter != 0:
                loopingcounter: int
                for loopingcounter in range(1, len(featurefileheaderlist)):
                    # row[0] is the id
                    # header is _id	confidence	count	Pathway	module	Gene	Transcript
                    #print("*****",row[0],"\t", featurefileheaderlist[loopingcounter],"\t", row[1:][loopingcounter])

                    featureinfodict2d[row[0]][featurefileheaderlist[loopingcounter]] = row[1:][loopingcounter]
            elif featurecounter == 0:
                # header row
                for sample in row[1:]:
                    featurefileheaderlist.append(sample)
                    featureinfodict2d[sample] = {}
            featurecounter = featurecounter + 1
    print("After reading in sample info records ", samplecounter)
    ##################################################################################################
    return sequencedict, featureinfodict2d, datainfodict2d, missingsamples, subjectstreatedwithdrugs


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
        sequencedict, featureinfodict2d, datainfodict2d, missingsamples, subjectstreatedwithdrugs = findgenesofsignificance(datafilename, featuresfilename, samplesfilename, sequencesfilename)
        print("After reading in samples, data features")
        for key in datainfodict2d.keys():
            print(key)
            print(datainfodict2d.get(key))
    finally:
        print(list(checkforvalidparentheses("{([((()])))[][[()]]}")))
        print("After printing the first 15 lines")
