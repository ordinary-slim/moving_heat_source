import re

integerKeys = [
        "nels",
        "maxIter",
        "plot",
        "timeIntegration",
        "sourceTerm",
        "isAdvection",
        ]

def parseInput(lines):
    dic = {}
    for l in lines:
        pair = l.split()
        dic[pair[0]] = float(pair[1])

    for iK in integerKeys:
        if iK in dic:
            dic[iK] = int(dic[iK])
    return dic

def formatInputFile( fileName ):
    lines = []
    with open( fileName, 'r') as f:
        lines = f.readlines()
    for l in lines:
        l.rstrip("\n")
    return lines
