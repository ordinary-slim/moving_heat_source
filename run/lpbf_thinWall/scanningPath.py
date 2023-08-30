import yaml

params = {}
with open("input.yaml", 'r') as paramsFile:
    params = yaml.safe_load(paramsFile)

partLen = params["partLen"]
layerThickness = params["layerThickness"]
speed = max(params["HeatSourceSpeed"])

def writeGcode(gcodeFile="Path.gcode", nLayers=None):
    gcodeLines = []
    # Start in -Y
    X = 0.0
    Y = - partLen[1] / 2
    Z = 0.0
    E = 0.0
    numLayers = int(partLen[2] / layerThickness)
    if nLayers is not None:
        numLayers = min( numLayers, nLayers )
    
    gcodeLines.append( "G0 F{}".format(speed, X, Y, Z) )
    for ilayer in range(numLayers):
        Z = layerThickness * (ilayer+1)
        E += 0.1
        gcodeLines.append( "G0 X{} Y{} Z{:.2f}".format(X, Y, Z) )
        Y = -Y
        gcodeLines.append( "G1 Y{} E{:.2f}".format(Y, E) )

    with open(gcodeFile, 'w') as gfile:
        gfile.writelines( [line+"\n" for line in gcodeLines] )

if __name__=="__main__":
    writeGcode("Path.gcode")
