import re
import argparse

def getTimes(pvdLines):
    times = []
    for idx, line in enumerate(pvdLines):
        match = re.search(r'timestep="(.*)" group', line)
        if match:
            times.append( (idx, float(match.group(1))) )
    return times

def main(rawPvdFile, targetPvdFile):
    with open(rawPvdFile, 'r') as rf:
        rawPvdLines = rf.readlines()
    with open(targetPvdFile, 'r') as tf:
        targetPvdLines = tf.readlines()
    rawTimes = getTimes(rawPvdLines)
    targetTimes = getTimes(targetPvdLines)

    indices2keep = []
    currentIdxRaw = 0
    for _, time in targetTimes:
        for idx, rtime in rawTimes[currentIdxRaw:]:
            if abs(time - rtime) < 1e-10:
                indices2keep.append(idx)
                currentIdxRaw = idx + 1
                break

    lines2write = [rawPvdLines[idx] for idx in indices2keep]
    nameRaw = rawPvdFile.split(".")[0]
    nameTarget = targetPvdFile.split(".")[0]
    nameFiltered = nameRaw + "_filtered2_" + nameTarget + ".pvd"
    lines2write = rawPvdLines[0:3] + lines2write + rawPvdLines[-2:]

    with open(nameFiltered, 'w') as nf:
        nf.writelines( lines2write )
    print("Wrote " + nameFiltered )


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_file')
    parser.add_argument('target_file')
    args = parser.parse_args()
    main(rawPvdFile=args.raw_file, targetPvdFile=args.target_file)
