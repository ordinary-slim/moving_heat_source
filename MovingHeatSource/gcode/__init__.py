import numpy as np
from enum import IntEnum
from MovingHeatSource.cpp import TrackType

def gcode2laserPath( gcodeFile, defaultPower=100.0 ):
    class Index(IntEnum):
        X = 0
        Y = 1
        Z = 2

    gf = open( gcodeFile, 'r' )

    coordinates = []
    times = []
    speeds = []
    powers = []
    trackTypes = []

    currentPosition = np.array([0.0, 0.0, 0.0])
    currentTime = 0.0
    currentSpeed = 10
    currentPower = 0.0

    for rawline in gf:
        hasMotion = False
        trackType = TrackType.cooling
        line = rawline.rstrip("\n")
        line = line.split(";", 1)[0]# remove comments

        instructions = line.split()
        for instruction in instructions:
            instructionType = instruction[0]
            instructionValue = float(instruction.lstrip(instructionType))
            if instructionType == "G":
                if instructionValue == 4:
                    trackType = TrackType.dwelling
            elif instructionType == "F":
                currentSpeed = instructionValue
            elif instructionType in ["X", "Y", "Z"]:
                currentPosition[ int( Index[instructionType] ) ] = instructionValue
                hasMotion = True
            elif instructionType == "E":
                if instructionValue > 0.0:
                    trackType = TrackType.printing
            elif instructionType == "P":
                if trackType in [TrackType.dwelling, TrackType.recoating]:
                    currentTime += instructionValue
            elif instructionType == "R":
                if (trackType == TrackType.dwelling) and (instructionValue > 0):
                    trackType = TrackType.recoating

        if hasMotion:
            if (len(coordinates) > 0 ):
                currentTime = times[-1] + np.linalg.norm(currentPosition - coordinates[-1]) / currentSpeed
        if trackType == TrackType.printing:
            currentPower = defaultPower
        else:
            currentPower = 0.0

        if hasMotion:
            coordinates.append( currentPosition.copy() )
            times.append( currentTime )
            speeds.append( currentSpeed )
            powers.append( currentPower )
            trackTypes.append( trackType )
        elif trackType in [TrackType.dwelling, TrackType.recoating]:
            coordinates.append( currentPosition.copy() )
            times.append( currentTime )
            speeds.append( 0.0 )
            powers.append( 0.0 )
            trackTypes.append( trackType )

    gf.close()

    return coordinates, times, speeds, powers, trackTypes
