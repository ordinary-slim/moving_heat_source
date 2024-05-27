import re
import pdb
import subprocess
import time

inputFile = "input.yaml"

def rewrite_input( elsPerRadius, tstepsPerRadius ):
    # Read file
    lines = []
    with open( inputFile, "r" ) as ipfile:
        lines = ipfile.readlines()
    # Replace
    for idx, line in enumerate(lines):
        if "elsPerRadius" in line:
            lines[idx] = re.sub(r"\d+", str(elsPerRadius), line)
        elif "tstepsPerRadius" in line:
            lines[idx] = re.sub(r"\d+", str(tstepsPerRadius), line)
    with open( inputFile, "w" ) as ipfile:
        ipfile.writelines(lines)

def write_runsh( coupledRun, elsPerRadius, tstepsPerRadius ):
    caseName = None
    runType = None
    if coupledRun:
        caseName = "coupled_elsPerRad{}".format(
                elsPerRadius )
        runType = "coupled"
    else:
        caseName = "reference_elsPerRad{}_tstepsPerRad{}".format(
                elsPerRadius, tstepsPerRadius )
        runType = "reference"

    # Write run.sh
    runSh = [
    "#!/bin/bash",
    "#SBATCH --job-name={}".format( caseName ),
    "#SBATCH --output=output-job_%j.out",
    "#SBATCH --error=output-job_%j.err",
    "#SBATCH --partition=HighParallelization",
    "#SBATCH --ntasks=1",
    "#SBATCH --mem=30G",
    "#SBATCH --time=4-0",
    "########### Further details -> man sbatch ##########",
    "source /data0/home/mslimani/moving_heat_source/modules.sh",
    "source /data0/home/mslimani/moving_heat_source/.venv/bin/activate",
    "python3 main.py --run-{} --case-name={}".format(
        runType,
        caseName,
        ),
    ]

    fileName = caseName + ".sh"
    with open(fileName, "w") as runsh:
        runsh.writelines( line + "\n" for line in runSh )
    return fileName

def checkSlurmJobState( jobid ) :
    checkCmd = "seff {}".format( jobid )
    checkOut = subprocess.check_output(checkCmd, shell=True)
    acceptableStates = ['RUNNING', 'PENDING', 'COMPLETED']
    for line in checkOut.decode().split("\n"):
        if line.startswith('State'):
            if any([s in line for s in acceptableStates]):
                print( "{} is one of the following: {}.".format( jobid, acceptableStates ) )
                return 0
    print( "Something went wrong with {}!".format( jobid ) )
    return 1

def getSbatchJobId( submitStr ):
    return float(re.search("job (\d+)", submitStr).group(1))

def runSbatch( sbatchFile ):
    runCmd = "sbatch {}".format( sbatchFile ) 
    submitStr = subprocess.check_output(runCmd, shell=True)
    jobid = getSbatchJobId( submitStr.decode() )
    time.sleep( 10 )
    checkSlurmJobState(jobid)
    return jobid

def main():
    # Loop over desired properties
    for elsPerRadius in [2, 4, 8, 16]:
        for tstepsPerRadius in [2, 4, 8, 16]:
            if tstepsPerRadius > elsPerRadius:
                continue
            rewrite_input( elsPerRadius, tstepsPerRadius )
            fileName = write_runsh( False, elsPerRadius, tstepsPerRadius )
            runSbatch( fileName )

    for elsPerRadius in [2, 4, 8, 16]:
        rewrite_input( elsPerRadius, 2 )
        fileName = write_runsh( True, elsPerRadius, 2 )
        runSbatch( fileName )

if __name__=="__main__":
    main()
