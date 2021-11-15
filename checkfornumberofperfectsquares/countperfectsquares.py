# Author Lalitha Viswanathan
# Code to find number of perfect squares between two integers min and max
import math
from functools import lru_cache
@lru_cache(maxsize=32)

def findperfectsquares(minnum: int, maxnum: int):
    """

    :type minnum: object
    """
    numberofperfectsquares: int = 0
    for num in range(minnum, maxnum+1):
        if checkifnumberisperfectsquare(num) is True:
            numberofperfectsquares = numberofperfectsquares + 1

    print("Number of perfect squares ", numberofperfectsquares)
    return numberofperfectsquares

# Author Lalitha Viswanathan
# To check if a given number is a perfect square
# square root and then square should give the number to be checked
def checkifnumberisperfectsquare(numbertobechecked: int) -> bool:
    if int((math.sqrt(numbertobechecked)) + 0.5) ** 2 == numbertobechecked:
        return True
    else:
         return False
    fib.cache_clear()