import numpy as np
from enum import IntEnum

def gcode2laserPath( gcodeFile, defaultPower=100.0 ):
    class Index(IntEnum):
        X = 0
        Y = 1
        Z = 2

    gf = open( gcodeFile, 'r' )

    coordinates = []
    speeds = []
    powers = []
    arePrinting = []

    currentPosition = np.array([0.0, 0.0, 0.0])
    currentSpeed = 10
    currentPower = 0.0

    for rawline in gf:
        hasCoordinate = False
        hasPrinting = False
        line = rawline.rstrip("\n")
        line = line.split(";", 1)[0]# remove comments

        instructions = line.split()
        for instruction in instructions:
            instructionType = instruction[0]
            instructionValue = float(instruction.lstrip(instructionType))
            if instructionType == "F":
                currentSpeed = instructionValue
            elif instructionType in ["X", "Y", "Z"]:
                currentPosition[ int( Index[instructionType] ) ] = instructionValue
                hasCoordinate = True
            elif instructionType == "E":
                if instructionValue > 0.0:
                    hasPrinting = True

        if hasPrinting:
            currentPower = defaultPower
        else:
            currentPower = 0.0

        if hasCoordinate:
            coordinates.append( currentPosition.copy() )
            speeds.append( currentSpeed )
            powers.append( currentPower )
            arePrinting.append( int(hasPrinting) )

    gf.close()

    return coordinates, speeds, powers, arePrinting
