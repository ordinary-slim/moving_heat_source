import MovingHeatSource as mhs
from MovingHeatSource.adaptiveStepper import LpbfAdaptiveStepper, deactivateBelowSurface

class CustomStepper(LpbfAdaptiveStepper):
    def setBCs( self ):
        self.pFixed.setConvection( resetBcs = False )
        if self.isCoupled:
            self.pMoving.setConvection( resetBcs = False )

def runCoupled(caseName="fixed"):
    from meshing import getMeshPhysical, fineElSize
    from scanningPath import writeGcode

    inputFile = "input.yaml"
    problemInput = mhs.readInput( inputFile )
    writeGcode( problemInput["path"] )
    # read input
    problemInput = mhs.readInput( inputFile )
    fixedProblemInput = dict( problemInput )

    # Mesh
    meshFixed = getMeshPhysical()

    pFixed = mhs.Problem(meshFixed, fixedProblemInput, caseName=caseName)

    deactivateBelowSurface( pFixed )

    driver = CustomStepper( pFixed,
                            maxAdimtDt=2.5,
                            elementSize=fineElSize,
                            threshold=0.3,
                            slowDown=True,
                            slowAdimDt=1,
                           )
    
    while not(driver.pFixed.mhs.path.isOver( driver.getTime() ) ) :
        driver.iterate()

def test():
    runCoupled(caseName="coupled")
    newFixed = "post_coupled/coupled_30.vtu"
    newMoving = "post_coupled_moving/coupled_moving_30.vtu"
    refFixed = "post_coupled_reference/coupled_30.vtu"
    refMoving = "post_coupled_moving_reference/coupled_moving_30.vtu"
    for ds1, ds2 in zip([newFixed, newMoving], [refFixed, refMoving]):
        assert mhs.meshio_comparison( ds1, ds2 )

if __name__=="__main__":
    test()
