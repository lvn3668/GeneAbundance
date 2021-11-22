# Features file contains metadata about the features whose abundance measures are in data.txt
# Metadata about features include confidence, count, pathways, modules, associated genes and transcript ids
# _id	confidence	count	Pathway	module	Gene	Transcript
import scipy


def checkforgaussiandistribution(data: list[float]) -> tuple[float, bool]:
    global checknormality
    # check for gaussian distribution in input dataset
    stat, p = scipy.stats.shapiro(data)
    #print("After calculating p value ", p)
    #print('stat=%.3f, p=%.3f' % (stat, p))
    isnormal: bool
    if p > 0.01:
            #print('Probably Gaussian')
            isnormal = True
    else:
            #print('Probably not Gaussian')
            isnormal = False
    return p, isnormal



