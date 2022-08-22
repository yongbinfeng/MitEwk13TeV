"""
script to generate and submit condor jobs
"""
import os
import time

def GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron=False, do5TeV=False):
    # creat the outdir first
    os.system('/usr/bin/eos root://cmseos.fnal.gov rm -r ' + outdir)
    os.system('/usr/bin/eos root://cmseos.fnal.gov mkdir -p ' + outdir)
    pwd = os.environ['PWD'];
    jobname = pwd + "/log/run_" + logsuffix + "_" + str(ijob)
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

        cmd = "root -l -q -b "
        cmd += macro
        cmd += '+\(\\"' + confname + '\\",\\"./' 
        if not do5TeV:
            if doElectron:
                # for electrons, do scale variations, with 1 sigma
                cmd += '\\",1,0,0,1,1,0,'+str(ijob) + '\)'
            else:
                cmd += '\\",0,0,1,1,0,'+str(ijob) + '\)'
        else:
            if doElectron:
                # 5TeV Electron channel
                cmd += '\\",1,0,0,0,1,0,'+str(ijob) + '\)'
            else:
                cmd += '\\",0,0,0,1,0,'+str(ijob) + '\)'
        bashscript.write(cmd+"\n")

        #bashscript.write("/usr/bin/eos root://cmseos.fnal.gov mkdir -p " + outdir + "\n")    
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
    time.sleep(0.1)



if __name__ == "__main__":
    ## Z -> mumu
    macro = "selectZmm.C"
    confname = "zmm_13.conf"
    outdir = "/store/user/yofeng/LowPU_Selection_pT20/Zmumu/"
    logsuffix = "selectZmm_zmm_13"
    if 0:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix)

    ### W -> munu
    macro = "selectWm.C"
    confname = "wm_13.conf"
    outdir = "/store/user/yofeng/LowPU_Selection/Wmunu"
    logsuffix = "selectWm_wm_13"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix)

    ## Z -> ee
    macro = "selectZee.C"
    confname = "zee_13.conf"
    outdir = "/store/user/yofeng/LowPU_Selection_pT20/Zee"
    logsuffix = "selectZee_zee_13"
    if 0:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)

    ## W -> enu
    macro = "selectWe.C"
    confname = "we_13.conf"
    outdir = "/store/user/yofeng/LowPU_Selection/Wenu"
    logsuffix = "selectWe_we_13"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)

    ## QCD control region (anti-isolated region) for W->enu
    macro = "selectAntiWe.C"
    confname = "we_13.conf"
    outdir = "/store/user/yofeng/LowPU_Selection/AntiWenu"
    logsuffix = "selectAntiWe_we_13"
    if 0:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)

    # 5 TeV
    # Z -> mumu
    macro = "selectZmm.C"
    confname = "zmm_5.conf"
    outdir = "/store/user/yofeng/LowPU_5TeV_Selection_pT20/Zmumu"
    logsuffix = "selectZmm_zmm_5"
    if 1:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, do5TeV = True)

    ## W -> munu
    macro = "selectWm.C"
    confname = "wm_5.conf"
    outdir = "/store/user/yofeng/LowPU_5TeV_Selection/Wmunu"
    logsuffix = "selectWm_wm_5"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, do5TeV = True)

    # Z -> ee
    macro = "selectZee.C"
    confname = "zee_5.conf"
    outdir = "/store/user/yofeng/LowPU_5TeV_Selection_pT20/Zee"
    logsuffix = "selectZee_zee_5"
    if 1:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True, do5TeV = True)

    # W -> enu
    macro = "selectWe.C"
    confname = "we_5.conf"
    outdir = "/store/user/yofeng/LowPU_5TeV_Selection/Wenu"
    logsuffix = "selectWe_we_5"
    if 1:
        for ijob in range(10):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True, do5TeV = True)

    ## QCD control region (anti-isolated region) for W->enu
    macro = "selectAntiWe.C"
    confname = "we_5.conf"
    outdir = "/store/user/yofeng/LowPU_5TeV_Selection/AntiWenu"
    logsuffix = "selectAntiWe_we_5"
    if 0:
        for ijob in range(10):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)
