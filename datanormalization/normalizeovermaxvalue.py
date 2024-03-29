# Author: Lalitha Viswanathan
# Package to normalize input data over maximum
from readdatafile.readdatafile import normalize

def normalize_over_maximum_value(dicttobenormalized: dict):
    ###################################################################################################
    # Normalize the featureid to count association
    maxcount = max(dicttobenormalized.values())
    dicttobenormalized = {k: normalize(v, maxcount) for k, v in
                                                             dicttobenormalized.items()}
    ###################################################################################################
    return dicttobenormalized