# Author: Lalitha Viswanathan
# Function to parse nested parenthesis


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
