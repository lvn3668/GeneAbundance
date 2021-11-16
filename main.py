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
import numpy
from sklearn.preprocessing import LabelEncoder

from associationbetweenfeatureconfidenceandcount.associationbetweenconfidenceandcount import \
    findassociationbetweentranscriptconfidenceandtranscriptcount
from associationbetweenlengthanddaysinpool.associationbetweenlengthanddaysinpool import \
    findassociationbetweensequencelengthanddaysinpool
from checkfornumberofperfectsquares.countperfectsquares import findperfectsquares
from findgenesofsignificancebasedonabundancemeasures.findgenesofsignificance import \
    readmetadataaboutsamplesfeaturesandexpressionvalues
from genesrankedbypathwayrank.generankingbasedonpathwayrank import rankgenesbasedonpathways
from normalizeexpressionvalues.normalizeexpressionmeasures import normalizeexpressionmeasures
from reversecomplementDNA.reversecomplementDNA import reversecomplementDNA, validationofreversecomplement
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import miceforest as mf
from sklearn.datasets import load_iris
import pandas as pd
import numpy as np


# Explore the `samples` metadata. Calculate how many sequenced samples are available
# for each subject (\"Subject\" column) across treatment usage (\"Treated_with_drug\" column).
# Hint: You can use a \"groupby\" function or a \"table\" function.**"

# Transform the `data` table to relative abundances - this is done by dividing the counts in each sample and
# each feature by the sum of the counts for all features in each sample. Name this new dataframe `data_ra`.**"
# Determine the top 10 genes with the highest mean relative abundance found across both subjects.
# The Gene label can be found in the `features` dataframe.**"

def findgenesofsignificance(datafilename: str, featuresfilename: str, samplesfilename: str, sequencesfilename: str):
    sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength, \
    sample_expressedfeature_or_transcript_expressionvalues_matrix, \
    meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
    expressedtranscriptspersample_sampletotranscriptassoc, sample_expressedfeature_or_transcript_expressionvalues_matrix, \
    samplid_to_subject_association, \
    subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
    sampleid_to_dayssinceexperimentstarted_normalized, features_metadata, transcript_or_feature_to_gene_association, \
    transcript_or_feature_to_confidence_association, \
    transcript_or_feature_to_count_association, gene_to_transcript_or_feature_association, \
    pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
    gene_to_pathway_association, pathway_to_gene_association, sample_to_gene_association, gene_to_sample_association, \
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues = readmetadataaboutsamplesfeaturesandexpressionvalues(
        datafilename, featuresfilename, samplesfilename, sequencesfilename)

    genesrankedbypathwayexpression = rankgenesbasedonpathways(geneid_expressionvalues,
                                                              gene_to_pathway_association,
                                                              pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes)
    import json
    print(json.dumps(genesrankedbypathwayexpression, sort_keys=True, indent=2))
    print("After printing gene list ranked by pathway expression NOT normalized ( num genes num modules correlated)")
    genesmissingpathwayassociations: list[str] = []
    for gene in geneid_expressionvalues.keys():
        # Multiply by confidence to normalize
        # Use counts to normalize
        geneid_expressionvalues[gene] = normalizeexpressionmeasures(gene, geneid_expressionvalues,
                                                                    genesrankedbypathwayexpression,
                                                                    gene_to_transcript_or_feature_association,
                                                                    transcript_or_feature_to_confidence_association,
                                                                    transcript_or_feature_to_count_association,
                                                                    gene_to_sample_association,
                                                                    sampleid_to_dayssinceexperimentstarted_normalized,
                                                                    sampleid_to_length_association_normalizedoverrmaxlength,
                                                                    genesmissingpathwayassociations)

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

        print("Gene ", gene, " Relative abundances ", geneid_expressionvalues.get(gene))

    # Determine the top 10 genes with the highest mean relative abundance found across both subjects.
    # The Gene label can be found in the `features` dataframe.**"
    ##################################################################################################

    return sequence_header_to_fnasequence_association, features_metadata, \
           sample_expressedfeature_or_transcript_expressionvalues_matrix, \
           samples_with_missing_expression_values, subjects_treatedwithdrugs_to_sample_association, \
           sampleid_to_length_association_normalizedoverrmaxlength, \
           sampleid_to_dayssinceexperimentstarted_normalized, geneid_expressionvalues, \
           transcript_or_feature_to_confidence_association, \
            transcript_or_feature_to_count_association




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
        sequencedict, featureinfodict2d, datainfodict2d, missingsamples, subjectstreatedwithdrugs, \
        sampleid_to_length_association_normalizedoverrmaxlength, sampleid_to_dayssinceexperimentstarted_normalized, \
        geneid_expressionvalues, transcript_or_feature_to_confidence_association, \
        transcript_or_feature_to_count_association = findgenesofsignificance(
                datafilename, featuresfilename, samplesfilename, sequencesfilename)
        # 8 samples have both days in pool and seq length
        # for those 8 samples,
        findassociationbetweensequencelengthanddaysinpool(sampleid_to_length_association_normalizedoverrmaxlength,
        sampleid_to_dayssinceexperimentstarted_normalized)
        findassociationbetweentranscriptconfidenceandtranscriptcount(transcript_or_feature_to_confidence_association,
        transcript_or_feature_to_count_association)
        print("Exit Before imputation")
        exit(1)
        from scipy.stats import beta
        #a, b = 1., 2.
        #x = beta.rvs(a, b, size=1000)
        #a1, b1, loc1, scale1 = beta.fit(x)
        #print("******* MLE fit ", a1, b1, loc1, scale1)
        # convert gene transcript expression vals to panda df
        df_genettranscript_exprvals = pd.DataFrame.from_dict(geneid_expressionvalues, orient='index').T
        print(df_genettranscript_exprvals)
        imp = SimpleImputer(missing_values=np.NaN, strategy='mean')
        imp.fit(df_genettranscript_exprvals)
        df_genettranscript_exprvals = imp.transform(df_genettranscript_exprvals)
        print(imp.transform(df_genettranscript_exprvals))
        np.random.seed(42)
        kf = KFold(n_splits=4)
        scores = []
        for train_index, test_index in kf.split(df_genettranscript_exprvals):
            X_train, X_test = df_genettranscript_exprvals[train_index], df_genettranscript_exprvals[test_index]
            y_train, y_test = df_genettranscript_exprvals[train_index], df_genettranscript_exprvals[test_index]

            clf = LinearRegression()
            clf.fit(X_train, y_train)
            y_test_pred = clf.predict(X_test)
            scores.append(r2_score(y_test, y_test_pred))

        print("Mice forest algorithm scores ", scores)

        # importing the dataset into kaggle
        df = pd.read_csv("data/test.csv")
        print(df.keys())

        # Load data and introduce missing values
        from sklearn import metrics
        from sklearn.model_selection import train_test_split

        le = LabelEncoder()
        df['Sex'] = le.fit_transform(df['Sex'])
        newdf = df

        y = df['SibSp']
        print(y)

        X_train, X_test, y_train, y_test = train_test_split(df, y, test_size=0.3)
        print(" training data ", X_train)
        print(" y axis training data ", y_train)
        from sklearn.linear_model import LogisticRegression

        lr = LogisticRegression()

        lr.fit(X_train, y_train)
        print("After lr.fit")
        pred = lr.predict(X_test)

        print(metrics.accuracy_score(pred, y_test))
        exit(1)

    finally:
        #print(list(checkforvalidparentheses("{([((()])))[][[()]]}")))
        print("After printing the first 15 lines")
