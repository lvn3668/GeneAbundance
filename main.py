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
# "Imagine a grocery store with multiple aisles and multiple people waiting in each line.
# You quickly look at the items in each customers basket and know exactly how long it will take them to check out.
# If each customer always chooses the shortest line, calculate how long it will take for all customers to go through checkout.
# For example if you have 3 customers who will take 5, 10, and 3 minutes each to go through checkout and there is only one lane available it will take 18 minutes,
# if there are two aisles available it will take ten minutes"
# Write a method to reverse translate the string “ACCTGGCCGTACCT”.
# Results should be “AGGTACGGCCAGGT”.
# Add a validation step to show that the results is what expected.\n",

from sklearn.preprocessing import LabelEncoder
from associationbetweenfeatureconfidenceandcount.associationbetweenconfidenceandcount import \
    findassociationbetweentranscriptconfidenceandtranscriptcount
from associationbetweenlengthanddaysinpool.associationbetweenlengthanddaysinpool import \
    findassociationbetweensequencelengthanddaysinpool
from checkfornumberofperfectsquares.countperfectsquares import findperfectsquares
from findgenesofsignificancebasedonabundancemeasures.findgenesofsignificance import \
    readmetadataaboutsamplesfeaturesandexpressionvalues
from genesrankedbypathwayrank.generankingbasedonpathwayrank import rankgenesbasedonpathways
from genexpressiondataimputation.imputegeneexpressiondata import imputegenexpressiondata
from normalizeexpressionvalues.normalizeexpressionmeasures import normalizeexpressionmeasures
from reversecomplementDNA.reversecomplementDNA import reversecomplementDNA, validationofreversecomplement
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import load_iris
import rpy2
from scipy.stats import beta
from sortedcontainers import SortedList


# Load data and introduce missing values

# Explore the `samples` metadata. Calculate how many sequenced samples are available
# for each subject (\"Subject\" column) across treatment usage (\"Treated_with_drug\" column).
# Hint: You can use a \"groupby\" function or a \"table\" function.**"

# Transform the `data` table to relative abundances - this is done by dividing the counts in each sample and
# each feature by the sum of the counts for all features in each sample. Name this new dataframe `data_ra`.**"
# Determine the top 10 genes with the highest mean relative abundance found across both subjects.
# The Gene label can be found in the `features` dataframe.**"

def findgenesofsignificance(datafilename: str, featuresfilename: str, samplesfilename: str, sequencesfilename: str):
    expressedtranscriptspersample_sampletotranscriptassoc: dict
    sequence_header_to_fnasequence_association, sampleid_to_length_association_normalizedoverrmaxlength, \
    sample_expressedfeature_or_transcript_expressionvalues_matrix, \
    meanexpression_per_sample_for_all_features_or_transcripts_sorted, \
    expressedtranscriptspersample_sampletotranscriptassoc, sample_expressedfeature_or_transcript_expressionvalues_matrix, \
    samplid_to_subject_association, \
    subjects_treatedwithdrugs_to_sample_association, samples_with_missing_expression_values, \
    sampleid_to_dayssinceexperimentstarted_normalized, features_metadata, transcript_or_feature_to_gene_association, \
    transcript_or_feature_to_confidence_normalized_association, \
    transcript_or_feature_to_count_normalized_association, gene_to_transcript_or_feature_association, \
    pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes, \
    gene_to_pathway_association, pathway_to_gene_association, sample_to_gene_association, gene_to_sample_association, \
    expressionvalues_subject_drugtreatment_sampleid_transcript_gene, geneid_expressionvalues, sample_treatedwithdrugs_association, \
    sample_transcript_expressionval_association, geneid_to_transcriptexprvals_with_featureids_association = readmetadataaboutsamplesfeaturesandexpressionvalues(
        datafilename, featuresfilename, samplesfilename, sequencesfilename)

    genesrankedbypathwayexpression = rankgenesbasedonpathways(geneid_expressionvalues,
                                                              gene_to_pathway_association,
                                                              pathways_ranked_by_positivecorrelnbetnnummodules_and_numgenes)

    genesmissingpathwayassociations: list[str] = []
    for gene in geneid_expressionvalues.keys():
        # Multiply by confidence to normalize
        # Use counts to normalize
        ####################################
        # Also return gene - transcript - measure
        geneid_expressionvalues[gene] = normalizeexpressionmeasures(gene, geneid_expressionvalues,
                                                                    genesrankedbypathwayexpression,
                                                                    gene_to_transcript_or_feature_association,
                                                                    transcript_or_feature_to_confidence_normalized_association,
                                                                    transcript_or_feature_to_count_normalized_association,
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
        # Pathway ranking (done)
        # build pathway module gene association
        # A gene can be implicated in multiple pathways and a pathway should have mutiple modules
        # A gene from a well expressed pathway -> HOW to normalize

        # to check: functional association and gene abundance expression

        # 11/19
        # Impute sequence length using correlation between length and num expressed transcripts / features
        # print("Gene ", gene, " Relative abundances ", geneid_expressionvalues.get(gene))

    # Determine the top 10 genes with the highest mean relative abundance found across both subjects.
    # The Gene label can be found in the `features` dataframe.**"
    ##################################################################################################

    ##############################################
    # gene expression values is 1: many
    # transcript to gene sanity check should be 1: 1
    # gene to trancript sanit check should be 1: many
    #
    combined_dict: dict = {}
    for gene in geneid_expressionvalues:
        if gene not in combined_dict.keys():
            combined_dict[gene] = {}
        combined_dict[gene]["NormalizedExprVals"] = geneid_expressionvalues.get(gene)
        combined_dict[gene]["NormalizedPathwayRank"] = genesrankedbypathwayexpression.get(gene)
        combined_dict[gene]["NormalizedSampleLength"] = []
        combined_dict[gene]["NormalizedDaysInPool"] = []
        for s in gene_to_sample_association.get(gene):
            combined_dict[gene]["NormalizedSampleLength"].append(
                sampleid_to_length_association_normalizedoverrmaxlength.get(s))
            combined_dict[gene]["NormalizedDaysInPool"].append(sampleid_to_dayssinceexperimentstarted_normalized.get(s))
            if expressedtranscriptspersample_sampletotranscriptassoc.get(s) is not None:
                combined_dict[gene]["NormalizedConfidence"] = []
                combined_dict[gene]["NormalizedCount"] = []
                combined_dict[gene]["TreatedWithDrugs"] = None
                for t in expressedtranscriptspersample_sampletotranscriptassoc.get(s):
                    combined_dict[gene]["NormalizedConfidence"].append(
                        transcript_or_feature_to_confidence_normalized_association.get(t))
                    combined_dict[gene]["NormalizedCount"].append(
                        transcript_or_feature_to_count_normalized_association.get(t))
                combined_dict[gene]["TreatedWithDrugs"] = sample_treatedwithdrugs_association.get(s)
    ##############################################
    # Dump of combined dict
    import json
    print(json.dumps(combined_dict, sort_keys=True, indent=2))

    return sequence_header_to_fnasequence_association, features_metadata, sample_expressedfeature_or_transcript_expressionvalues_matrix, samples_with_missing_expression_values, subjects_treatedwithdrugs_to_sample_association, sampleid_to_length_association_normalizedoverrmaxlength, sampleid_to_dayssinceexperimentstarted_normalized, geneid_expressionvalues, transcript_or_feature_to_confidence_normalized_association, transcript_or_feature_to_count_normalized_association


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
                                                          sampleid_to_dayssinceexperimentstarted_normalized)  # positive correlation
        findassociationbetweentranscriptconfidenceandtranscriptcount(transcript_or_feature_to_confidence_association,
                                                                     transcript_or_feature_to_count_association)

        # convert gene transcript expression vals to panda df
        df_genetranscript_exprvals = imputegenexpressiondata(geneid_expressionvalues)
        np.random.seed(42)
        # kf = KFold(n_splits=8)
        scores = []
        # for train_index, test_index in kf.split(df_genetranscript_exprvals):

        # X_train, X_test = df_genetranscript_exprvals[train_index], df_genetranscript_exprvals[test_index]
        # y_train, y_test = df_genetranscript_exprvals[train_index], df_genetranscript_exprvals[test_index]
        # clf = LinearRegression()
        # clf.fit(X_train, y_train)
        # y_test_pred = clf.predict(X_test)
        # scores.append(r2_score(y_test, y_test_pred))

        # test_size = 0.3
        # X_train, X_test, y_train, y_test = train_test_split(df_genetranscript_exprvals[train_index],
        # df_genetranscript_exprvals[test_index])
        # print(" training data ", X_train)
        # print(" y axis training data ", y_train)
        # lr = LogisticRegression()
        # lr.fit(X_train, y_train)
        # print("After lr.fit")
        # pred = lr.predict(X_test)
        # print(metrics.accuracy_score(pred, y_test))

        print("Mice forest algorithm scores ", scores)
        # le = LabelEncoder()

        import rpy2.robjects as robjects

        r = robjects.r

        m = r.matrix(r.rnorm(100), ncol=5)
        pca = r.princomp(m)
        r.plot(pca, main="Eigen values")
        r.biplot(pca, main="biplot")

        exit(1)

    finally:
        # print(list(checkforvalidparentheses("{([((()])))[][[()]]}")))
        print("After printing the first 15 lines")
