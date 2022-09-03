"""
script to generate and submit condor jobs
"""
import os
import time
from collections import OrderedDict

def GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs, doElectron=False, do5TeV=False):
    # creat the outdir first
    os.system('/usr/bin/eos root://cmseos.fnal.gov rm -r ' + outdir)
    os.system('/usr/bin/eos root://cmseos.fnal.gov mkdir -p ' + outdir)
    pwd = os.environ['PWD'];
    for ijob in xrange(njobs):
        jobname = pwd + "/log/run_" + logsuffix + "_samp" + str(isamp) + "_" + str(njobs) + "jobs_" + str(ijob) + "th"
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
            bashscript.write("cd $CMSSW_BASE/src/MitEwk13TeV/Selection\n")
            bashscript.write("\n")

            cmd_end = str(njobs) + ',' + str(ijob)+ ',' + str(isamp) + '\)'

            cmd = "root -l -q -b "
            cmd += macro
            cmd += '++\(\\"' + confname + '\\",\\"./' 
            if not do5TeV:
                if doElectron:
                    # for electrons, do scale variations, with 1 sigma
                    cmd += '\\",1,0,0,1,'
                else:
                    cmd += '\\",0,0,1,'
            else:
                if doElectron:
                    # 5TeV Electron channel
                    cmd += '\\",1,0,0,0,'
                else:
                    cmd += '\\",0,0,0,'
            cmd += cmd_end
            bashscript.write(cmd+"\n")

            bashscript.write("xrdcp -f ./ntuples/*root root://cmseos.fnal.gov/" + outdir + "/\n")

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

        os.system("condor_submit " + jobname + ".condor")
        time.sleep(0.01)


if __name__ == "__main__":
    ## Z -> mumu
    macro = "selectZmm.C"
    confname = "zmm_13.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zmumu_pT20/"
    logsuffix = "selectZmm_zmm_13"
    njobs = [20, 10, 3, 10, 3, 3, 3, 40, 3, 3, 10, 20]
    if 0:
        for isamp in range(12):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp])

    ### W -> munu
    macro = "selectWm.C"
    confname = "wm_13.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Wmunu/"
    logsuffix = "selectWm_wm_13"
    njobs = [20, 3, 3, 10, 3, 3, 3, 80, 3, 3, 10, 100, 20, 10]
    if 1:
        for isamp in range(14):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp])

    ## Z -> ee
    macro = "selectZee.C"
    confname = "zee_13.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Zee_pT20/"
    logsuffix = "selectZee_zee_13"
    njobs = [20, 10, 3, 10, 3, 3, 3, 40, 3, 3, 10, 20]
    if 0:
        for isamp in range(12):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], doElectron = True)

    ## W -> enu
    macro = "selectWe.C"
    confname = "we_13.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/Wenu/"
    logsuffix = "selectWe_we_13"
    njobs = [20, 3, 3, 10, 3, 3, 3, 80, 3, 3, 10, 100, 20, 10]
    if 1:
        for isamp in range(14):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], doElectron = True)

    ### QCD control region (anti-isolated region) for W->munu
    macro = "selectAntiWm.C"
    confname = "wm_13.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/AntiWmunu/"
    logsuffix = "selectAntiWm_wm_13"
    njobs = [20, 3, 3, 10, 3, 3, 3, 80, 3, 3, 10, 100, 20, 10]
    if 1:
        for isamp in range(14):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp])

    ## QCD control region (anti-isolated region) for W->enu
    macro = "selectAntiWe.C"
    confname = "we_13.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/13TeV/Selections/AntiWenu/"
    logsuffix = "selectAntiWe_we_13"
    njobs = [20, 3, 3, 10, 3, 3, 3, 80, 3, 3, 10, 100, 20, 10]
    if 1:
        for isamp in range(14):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], doElectron = True)

    # 5 TeV
    # Z -> mumu
    macro = "selectZmm.C"
    confname = "zmm_5.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Zmumu_pT20/"
    logsuffix = "selectZmm_zmm_5"
    njobs = [20, 10, 3, 3, 3, 3, 10, 10, 20]
    if 0:
        for isamp in range(9):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], do5TeV = True)

    ## W -> munu
    macro = "selectWm.C"
    confname = "wm_5.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Wmunu/"
    logsuffix = "selectWm_wm_5"
    njobs = [20, 10, 3, 3, 3, 3, 10, 10, 20]
    if 1:
        for isamp in range(9):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], do5TeV = True)

    # Z -> ee
    macro = "selectZee.C"
    confname = "zee_5.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Zee_pT20/"
    logsuffix = "selectZee_zee_5"
    njobs = [20, 10, 3, 3, 3, 3, 10, 10, 20]
    if 0:
        for isamp in range(9):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], doElectron = True, do5TeV = True)

    # W -> enu
    macro = "selectWe.C"
    confname = "we_5.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/Wenu/"
    logsuffix = "selectWe_we_5"
    njobs = [20, 10, 3, 3, 3, 3, 10, 10, 20]
    if 1:
        for isamp in range(9):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], doElectron = True, do5TeV = True)

    ### QCD control region (anti-isolated region) for W->munu
    macro = "selectAntiWm.C"
    confname = "wm_5.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/AntiWmunu/"
    logsuffix = "selectAntiWm_wm_5"
    njobs = [20, 10, 3, 3, 3, 3, 10, 10, 20]
    if 1:
        for isamp in range(9):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], do5TeV = True)

    ## QCD control region (anti-isolated region) for W->enu
    macro = "selectAntiWe.C"
    confname = "we_5.conf"
    outdir = "/store/user/yofeng/Ntuples_LowPU/5TeV/Selections/AntiWenu/"
    logsuffix = "selectAntiWe_we_5"
    njobs = [20, 10, 3, 3, 3, 3, 10, 10, 20]
    if 1:
        for isamp in range(9):
            GenerateExecutable(macro, confname, outdir, isamp, logsuffix, njobs[isamp], doElectron = True, do5TeV = True)
