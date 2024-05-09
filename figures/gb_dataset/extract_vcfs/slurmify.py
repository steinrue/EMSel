import sys
import os
import pathlib2



TMP_DIR = 'slurmify_tmp'
NUM_CORES = 1
NUM_GB = 16




def readCommands (commandFile):

    # read file
    ifs = open (commandFile, 'r')
    theLines = ifs.readlines()
    ifs.close()

    # go through and clean
    cleaned = []
    for thisLine in theLines:
        cleaned.append (thisLine.strip())

    return cleaned


def prepareSlurmScripts (commandList, metaSubmitScript, tmpDir):

    # set up the tmp dir
    os.makedirs (tmpDir)

    # set up meta script
    metafs = open (metaSubmitScript, 'w')

    # go through commands
    for (cmdIdx, thisCmd) in enumerate (commandList):

        # generic name
        thisName = f"slurm{cmdIdx}"
        thisScriptName = pathlib2.Path ( f"{tmpDir}/{thisName}.sbatch")
        ofs = open (thisScriptName, "w")

        # needs to be in script
        ofs.write ("#!/bin/bash\n")
        # give it a name
        ofs.write (f"#SBATCH --job-name={thisName}\n")
        # only 1 node
        ofs.write ("#SBATCH --nodes=1\n")
        # only 1 task
        ofs.write ("#SBATCH --ntasks=1\n")
        # we want this many cpus
        ofs.write (f"#SBATCH --cpus-per-task={NUM_CORES}\n")
        # also a bit of time
        ofs.write ("#SBATCH --time=12:00:00\n")
        # some memory
        ofs.write (f"#SBATCH --mem={NUM_GB}gb\n")
        # output
        ofs.write (f"#SBATCH --output={tmpDir}/%x.o%j\n")
        ofs.write (f"#SBATCH --error={tmpDir}/%x.e%j\n")
        # and finally the command
        ofs.write (thisCmd + "\n")

        # be nice
        ofs.close()

        # add something to a meta script
        metafs.write ("sbatch %s\n" % thisScriptName)
    
    # close and then we done?
    metafs.close()


def main():

    # need exactly one additional argument
    assert (len(sys.argv) == 2)

    # get the commands
    commandList = readCommands (sys.argv[-1])

    prepareSlurmScripts (commandList, 'metasubmit.sh', TMP_DIR)


if __name__ == "__main__":
    main()
