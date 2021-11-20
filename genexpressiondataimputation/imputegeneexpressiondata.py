# Author Lalitha Viswanathan
# Package to impute gene expression data
import numpy as np
import pandas as pd
from pandas import DataFrame
from sklearn.impute import SimpleImputer


def imputegenexpressiondata(geneid_expressionvalues: dict) -> DataFrame:
    df_genetranscript_exprvals: DataFrame = pd.DataFrame.from_dict(geneid_expressionvalues, orient='index').T
    print(df_genetranscript_exprvals)
    imp = SimpleImputer(missing_values=np.NaN, strategy='mean')
    imp.fit(df_genetranscript_exprvals)
    df_genetranscript_exprvals = imp.transform(df_genetranscript_exprvals)
    return df_genetranscript_exprvals