def reversecomplementDNA(DNAstring: str) -> str:
    """
    :rtype: object
    :return: 
    :param DNAstring:
    :return:
    """
    print("DNA String ", DNAstring)
    if len(DNAstring) == 0:
        raise Exception("Invalid DNA String")
    stringarr: list[str] = list(DNAstring)
    reversecomplement: str = ""
    for char in stringarr:
        if char.upper() == 'A':
            reversecomplement = reversecomplement + "T"
        elif char.upper() == 'T':
            reversecomplement = reversecomplement + "A"
        elif char.upper() == 'G':
            reversecomplement = reversecomplement + "C"
        else:
            reversecomplement = reversecomplement + "G"

    reversecomplement = reversecomplement[::-1]
    return reversecomplement


def validationofreversecomplement(DNAstring: str, reversecomplementseq: str) -> bool:
    """
    :return:
    :param DNAstring:
    :param reversecomplementseq:
    :return:
    """
    if DNAstring == reversecomplementDNA(reversecomplementseq):
        return True
    else:
        return False