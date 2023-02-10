"""
script to generate and submit condor jobs
"""
import os
import time

def GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron=False, do5TeV=False):
    pwd = os.environ['PWD'];
    jobname = pwd + "/log/run_" + logsuffix + "_" + str(ijob)
    with open(jobname + ".sh", "w") as bashscript:
        bashscript.write("#!/bin/bash\n")
        bashscript.write("\n")
        bashscript.write("cd /afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Selection\n")
        bashscript.write("eval $(scramv1 runtime -sh);\n")
        bashscript.write("\n")

        # make the temp directory, copy the macro and libraries there
        bashscript.write("# create temp directory and copy macro and config there\n")
        tmpname = "/tmp/" + logsuffix + "/" + str(ijob) + "/"
        bashscript.write("mkdir -p " + tmpname + "\n")
        bashscript.write("cp rootlogon.C " + tmpname + "\n")
        bashscript.write("cp " + macro + " " + tmpname + "\n")
        bashscript.write("cp " + confname + " " + tmpname + "\n")
        bashscript.write("cd " + tmpname + "\n")
        bashscript.write("\n")

        cmd = "root -l -q "
        cmd += macro
        cmd += '++O\(\\"' + confname + '\\",\\"' + outdir
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
        bashscript.write(cmd)

    # change it to executable
    os.system("chmod +x " + jobname + ".sh")

    job_desc = """Universe = vanilla
Executable = {jobname}.sh
Log        = {jobname}.log
Output     = {jobname}.out
Error      = {jobname}.error
getenv      = True
environment = "LS_SUBCWD={here}"
+JobFlavour = "workday"
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
    outdir = "/eos/user/y/yofeng/LowPU/Selection_pT20/Zmumu/"
    logsuffix = "selectZmm_zmm_13"
    if 1:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix)

    ### W -> munu
    macro = "selectWm.C"
    confname = "wm_13.conf"
    outdir = "/eos/user/y/yofeng/LowPU/Selection/Wmunu"
    logsuffix = "selectWm_wm_13"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix)

    ## Z -> ee
    macro = "selectZee.C"
    confname = "zee_13.conf"
    outdir = "/eos/user/y/yofeng/LowPU/Selection_pT20/Zee"
    logsuffix = "selectZee_zee_13"
    if 1:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)

    ## W -> enu
    macro = "selectWe.C"
    confname = "we_13.conf"
    outdir = "/eos/user/y/yofeng/LowPU/Selection/Wenu"
    logsuffix = "selectWe_we_13"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)

    ## QCD control region (anti-isolated region) for W->enu
    macro = "selectAntiWe.C"
    confname = "we_13.conf"
    outdir = "/eos/user/y/yofeng/LowPU/Selection/AntiWenu"
    logsuffix = "selectAntiWe_we_13"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True)

    # 5 TeV
    # Z -> mumu
    macro = "selectZmm.C"
    confname = "zmm_5.conf"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/Selection_pT20/Zmumu"
    logsuffix = "selectZmm_zmm_5"
    if 1:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, do5TeV = True)

    ## W -> munu
    macro = "selectWm.C"
    confname = "wm_5.conf"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/Selection/Wmunu"
    logsuffix = "selectWm_wm_5"
    if 1:
        for ijob in range(14):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, do5TeV = True)

    # Z -> ee
    macro = "selectZee.C"
    confname = "zee_5.conf"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/Selection_pT20/Zee"
    logsuffix = "selectZee_zee_5"
    if 1:
        for ijob in range(12):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True, do5TeV = True)

    # W -> enu
    macro = "selectWe.C"
    confname = "we_5.conf"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/Selection/Wenu"
    logsuffix = "selectWe_we_5"
    if 1:
        for ijob in range(10):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True, do5TeV = True)

    ## QCD control region (anti-isolated region) for W->enu
    macro = "selectAntiWe.C"
    confname = "we_5.conf"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/Selection/AntiWenu"
    logsuffix = "selectAntiWe_we_5"
    if 1:
        for ijob in range(10):
            GenerateExecutable(macro, confname, outdir, ijob, logsuffix, doElectron = True, do5TeV = True)
