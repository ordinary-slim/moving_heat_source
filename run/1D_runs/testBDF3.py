import numpy as np
from prepos import *

def myBDF3():
    fileName = "testBDF3.txt"
    d = formatInputFile( fileName )
    d = parseInput( d )
    d["timeIntegration"] = 3
    d["nels"] = 5
    p = mhs.Problem()
    p.initialize( d )

    # write postfiles
    postFiles = []
    p.iterate()
    p.solution = 26*np.ones( p.solution.size )
    postFiles.append( writePost( p ) )
    p.iterate()
    p.solution = 27*np.ones( p.solution.size )
    postFiles.append( writePost( p ) )
    p.iterate()
    p.solution = 28*np.ones( p.solution.size )
    postFiles.append( writePost( p ) )

    # read postfiles
    initializeTimeIntegrator( p, postFiles )
    p.iterate()
    cppSol = p.solution
    print( "c++ sol= {}".format(cppSol) + "\n" )

    # try to match from python side
    dt = 1.0
    T = 28*np.ones( (6, 1) )
    Tn1 = 27*np.ones( (6, 1) )
    Tn2 = 26*np.ones( (6, 1) )
    M = np.diag( np.array([10, 20, 20, 20, 20, 10]) )
    K = np.array(
            [[0.5,-0.5,0,0,0,0],
            [-0.5,1,-0.5,0,0,0],
            [0,-0.5,1,-0.5,0,0],
            [0,0,-0.5,1,-0.5,0],
            [0,0,0,-0.5,1,-0.5],
            [0,0,0,0,-0.5,0.5]]
            )
    lhs = 0
    lhs += 11.0/6.0 * M / dt
    lhs += K
    print( "python lhs= {}".format(lhs) + "\n" )
    rhs = 0
    rhs += (+3.0   * M @ T) / dt
    rhs += (-1.5   * M @ Tn1) / dt
    rhs += (1.0/3.0* M @ Tn2) / dt
    pythonSol = np.linalg.solve( lhs, rhs )
    print( "python sol= {}".format(pythonSol) + "\n" )


if __name__=="__main__":
    myBDF3()

