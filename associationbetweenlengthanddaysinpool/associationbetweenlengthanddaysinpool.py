# Author: Lalitha Viswanathan
# Function to find association between length of fasta sequence associated with each sample and
# days in experimental pool
import numpy
import pandas as pd
from scipy.stats import stats

from checkforgaussiandistribution.gaussiandistributionusingshapirowilkis import checkforgaussiandistribution


def findassociationbetweensequencelengthanddaysinpool(sampleid_to_length_association_normalizedoverrmaxlength: dict,
                                                      sampleid_to_dayssinceexperimentstarted_normalized: dict,
                                                      ) -> object:
    df1 = pd.Series(data=sampleid_to_length_association_normalizedoverrmaxlength
                    ).to_frame().T
    df2 = pd.Series(sampleid_to_dayssinceexperimentstarted_normalized).to_frame().T
    indexcols = df1.columns.intersection(df2.columns)
    print("Index columns ", indexcols)
    intersected_df = pd.merge(df1, df2, on=list(df1.columns.intersection(df2.columns)))
    print("Intersected dataframe ", intersected_df)
    sample_length_daysinpool_association: dict[str, list] = {}
    for ky in list(
            sampleid_to_length_association_normalizedoverrmaxlength.keys() & sampleid_to_dayssinceexperimentstarted_normalized.keys()):
        sample_length_daysinpool_association[ky] = [sampleid_to_length_association_normalizedoverrmaxlength.get(ky),
                                                    sampleid_to_dayssinceexperimentstarted_normalized.get(ky)]
    df = pd.DataFrame.from_dict(sample_length_daysinpool_association).T

    df.columns = ["SequenceLength", "DaysInPool"]
    print(df)

    (pval_seqlength, normality_seqlength) = checkforgaussiandistribution(df.loc[:, "SequenceLength"])
    (pval_daysinpool, normality_daysinpool) = checkforgaussiandistribution(df.loc[:, "DaysInPool"])

    # Check for linear regression coefficients
    slope, intercept, r_value, p_value, std_err = stats.linregress(df.loc[:, "SequenceLength"], df.loc[:, "DaysInPool"])
    print(" Slope ", slope, " ", intercept, " ", r_value, " ", p_value, " ", std_err)

    # Calculate pearson's coefficient
    if normality_seqlength is True and normality_daysinpool is True:
        dataforfindingcovariance = list(
            zip(df.loc[:, "SequenceLength"], df.loc[:, "DaysInPool"]))
        print(dataforfindingcovariance)
        covariance_between_seqlength_to_daysinpool = numpy.cov(dataforfindingcovariance)
        print(covariance_between_seqlength_to_daysinpool)

    return covariance_between_seqlength_to_daysinpool