# Author Lalitha Viswanathan
# Function to find association between feature / transcript confidence and count
import pandas as pd
from scipy.stats import stats
from statsmodels.compat import numpy

from checkforgaussiandistribution.gaussiandistributionusingshapirowilkis import checkforgaussiandistribution


def findassociationbetweentranscriptconfidenceandtranscriptcount(
transcript_or_feature_to_confidence_association: dict,
        transcript_or_feature_to_count_association: dict):
    df1 = pd.Series(data=transcript_or_feature_to_confidence_association
                    ).to_frame().T
    df2 = pd.Series(transcript_or_feature_to_count_association).to_frame().T
    indexcols = df1.columns.intersection(df2.columns)
    print("Index columns ", indexcols)
    intersected_df = pd.merge(df1, df2, on=list(df1.columns.intersection(df2.columns)))
    print("Intersected dataframe ", intersected_df)
    transcript_confidence_count_association: dict[str, list] = {}
    for ky in list(
            transcript_or_feature_to_confidence_association.keys() &
            transcript_or_feature_to_count_association.keys()):
        transcript_confidence_count_association[ky] = [transcript_or_feature_to_confidence_association.get(ky),
                                                    transcript_or_feature_to_count_association.get(ky)]
    df = pd.DataFrame.from_dict(transcript_confidence_count_association).T

    df.columns = ["Confidence", "Count"]
    print(df)

    (pval_confidence, normality_confidence) = checkforgaussiandistribution(df.loc[:, "Confidence"])
    (pval_count, normality_count) = checkforgaussiandistribution(df.loc[:, "Count"])

    # Check for linear regression coefficients
    slope, intercept, r_value, p_value, std_err = stats.linregress(df.loc[:, "Confidence"], df.loc[:, "Count"])
    print(" Slope ", slope, " ", intercept, " ", r_value, " ", p_value, " ", std_err)

    # Calculate pearson's coefficient
    # If data is not normally distributed
    if normality_confidence is True and normality_count is True:
        dataforfindingcovariance = list(
            zip(df.loc[:, "Confidence"], df.loc[:, "Count"]))
        print(dataforfindingcovariance)
        covariance_between_seqlength_to_daysinpool = numpy.cov(dataforfindingcovariance)
        print(covariance_between_seqlength_to_daysinpool)
    else:
        print("Pvalue confidence normality confidence", pval_confidence, normality_confidence)
        print("Pvalue count normality count", pval_count, normality_count)
        covariance_between_seqlength_to_daysinpool = None


    return covariance_between_seqlength_to_daysinpool