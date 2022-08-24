"""
script to generate and submit condor jobs
"""
import os
import time
from collections import OrderedDict
import subprocess

def GenerateExecutable(macro, indir, outdir, logsuffix, is5TeV = False):

    os.system('/usr/bin/eos root://cmseos.fnal.gov rm -r ' + outdir)
    os.system('/usr/bin/eos root://cmseos.fnal.gov mkdir -p ' + outdir)

    # get the list of files in the indir
    cmd = '/usr/bin/eos root://cmseos.fnal.gov ls '
    cmd += indir
    infiles = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].split('\n')[:-1]

    pwd = os.environ['PWD'];
    if not is5TeV:
        dirname = pwd + "/log_13TeV/"
    else:
        dirname = pwd + "/log_5TeV/"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    idx = 0
    while idx < len(infiles):
        nfiles_per_job = 15
        if idx % nfiles_per_job == 0:
            ijob = idx / nfiles_per_job
            jobname = dirname + "run_" + logsuffix + "_" + str(ijob) + "job"
            with open(jobname + ".sh", "w") as bashscript:
                bashscript.write("#!/bin/bash\n")
                bashscript.write('\n')
                bashscript.write('echo "Starting job on " `date` #Date/time of start of job\n')
                bashscript.write('echo "Running on: `uname -a`" #Condor job is running on this node\n')
                bashscript.write('echo "System software: `cat /etc/redhat-release`" #Operating System on that node\n')
                bashscript.write('xrdcp -s root://cmseos.fnal.gov//store/user/yofeng/CMSSW_9_4_19.tgz .\n')
                bashscript.write('echo "copied CMSSW_9_4_19 from eos to local"\n')
                bashscript.write('\n')
                bashscript.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
                bashscript.write('tar -xf CMSSW_9_4_19.tgz\n')
                bashscript.write('rm CMSSW_9_4_19.tgz\n')
                bashscript.write('cd CMSSW_9_4_19/src/\n')
                bashscript.write('scramv1 b ProjectRename\n')
                bashscript.write('eval $(scram runtime -sh);\n')
                bashscript.write('scram b clean\n')
                bashscript.write('scram b -j 10\n')
                bashscript.write('echo "compiled"\n')
                bashscript.write('echo $CMSSW_BASE "is the CMSSW we have on the local worker node"\n')
                bashscript.write("cd $CMSSW_BASE/src/MitEwk13TeV/NtupleMod\n")
                bashscript.write("\n")

                while True:
                    fname = infiles[idx]
                    cmd = "root -l -q -b "
                    cmd += macro
                    cmd += '+\(\\"./ntuples/\\",\\"root://cmseos.fnal.gov/' + indir
                    if not is5TeV:
                        cmd += '\\",\\"13TeV\\",\\"' + fname + '\\"'
                    else:
                        cmd += '\\",\\"5TeV\\",\\"' + fname + '\\"'
                    cmd += '\)\n'
                    bashscript.write(cmd)

                    idx += 1
                    if idx % nfiles_per_job == 0 or idx >= len(infiles):
                        break

                bashscript.write("\nxrdcp -f ./ntuples/*root root://cmseos.fnal.gov/" + outdir + "/\n")

        # change it to executable
        os.system("chmod +x " + jobname + ".sh")

        job_desc = """Universe = vanilla
Executable = {jobname}.sh
Log        = {jobname}.log
Output     = {jobname}.out
Error      = {jobname}.error
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
queue 1\n
""".format(jobname = jobname, here = pwd)

        with open(jobname + ".condor", 'w') as outfile:
            outfile.write(job_desc)
            outfile.close()

        #os.system("condor_submit " + jobname + ".condor")
        #time.sleep(0.05)


if __name__ == "__main__":
    # W -> munu
    macro = "muonNtupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Wmunu/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/NtupleMod/Wmunu/"
    logsuffix = "muonNtupleMod_wm_13"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix)
    
    ## W -> enu
    macro = "eleNtupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Wenu/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/NtupleMod/Wenu/"
    logsuffix = "eleNtupleMod_we_13"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix)

    macro = "ZmmNTupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zmumu_pT20/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/NtupleMod/Zmumu/"
    logsuffix = "ZmmNTupleMod_zmumu_13"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix)

    macro = "ZeeNTupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zee_pT20/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/NtupleMod/Zee/"
    logsuffix = "ZeeNTupleMod_zee_13"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix)

    #
    # 5TeV results
    #
    # W -> munu
    macro = "muonNtupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Wmunu/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/NtupleMod/Wmunu/"
    logsuffix = "muonNtupleMod_wm_5"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix, True)

    ## W -> enu
    macro = "eleNtupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Wenu/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/NtupleMod/Wenu/"
    logsuffix = "eleNtupleMod_we_5"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix, True)

    macro = "ZmmNTupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Zmumu_pT20/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/NtupleMod/Zmumu/"
    logsuffix = "ZmmNTupleMod_zmumu_5"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix, True)

    macro = "ZeeNTupleMod.C"
    indir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Zee_pT20/"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/NtupleMod/Zee/"
    logsuffix = "ZeeNTupleMod_zee_5"
    if 1:
        GenerateExecutable(macro, indir, outdir, logsuffix, True)
