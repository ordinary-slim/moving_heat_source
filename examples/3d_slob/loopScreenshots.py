import subprocess

if __name__=="__main__":
    coupledDataSets = [
    "coupled_elsPerRad2.pvd",
    "coupled_elsPerRad4.pvd",
    "coupled_elsPerRad8.pvd",
    "coupled_elsPerRad16.pvd",]

    refDataSets = []
    for elsPerRad in [2, 4, 8, 16]:
        for tstepsPerRad in [2, 4, 8, 16]:
            if tstepsPerRad > elsPerRad:
                continue
            refDataSets.append(
                    "reference_elsPerRad{}_tstepsPerRad{}.pvd".format( elsPerRad, tstepsPerRad)
                    )

    for dataSet in coupledDataSets + refDataSets:
        subprocess.run(["pvpython", "screenshotContours.py", dataSet]) 
        subprocess.run(["pvpython", "screenshotError.py", dataSet]) 
    subprocess.run(["pvpython", "screenshotContours.py", "-i", "coupled_elsPerRad8_closeInterface.pvd"]) 
