"""
script to generate and submit condor jobs
"""
import os
import time

def GenerateExecutable(macro, indir, outdir, fname, logsuffix):
    pwd = os.environ['PWD'];
    jobname = pwd + "/log/run_" + logsuffix + "_" + str(fname)
    with open(jobname + ".sh", "w") as bashscript:
        bashscript.write("#!/bin/bash\n")
        bashscript.write("\n")
        bashscript.write("cd /afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV/NtupleMod\n")
        bashscript.write("eval $(scramv1 runtime -sh);\n")
        bashscript.write("\n")

        # make the temp directory, copy the macro and libraries there
        bashscript.write("# create temp directory and copy macro and config there\n")
        tmpname = "/tmp/" + logsuffix + "/" + str(fname) + "/"
        bashscript.write("mkdir -p " + tmpname + "\n")
        bashscript.write("cp rootlogon.C " + tmpname + "\n")
        bashscript.write("cp " + macro + " " + tmpname + "\n")
        bashscript.write("cd " + tmpname + "\n")
        bashscript.write("\n")

        cmd = "root -l -q "
        cmd += macro
        cmd += '++O\(\\"' + outdir + '\\",\\"' + indir
        cmd += '\\",\\"13\\",\\"' + fname + '.root\\"\)'
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
+JobFlavour = "tomorrow"
queue 1\n
""".format(jobname = jobname, here = pwd)

    with open(jobname + ".condor", 'w') as outfile:
        outfile.write(job_desc)
        outfile.close()

    os.system("condor_submit " + jobname + ".condor")
    time.sleep(0.1)


if __name__ == "__main__":
    # W -> munu
    macro = "muonNtupleMod.C"
    fnames = ["data_select.root", 
            "wm0_select.raw.root", "wm1_select.raw.root", "wm2_select.raw.root", 
            "wx0_select.raw.root", "wx1_select.raw.root", "wx2_select.raw.root", 
            "zxx_select.raw.root",
            "ww_select.raw.root", "wz_select.raw.root", "zz_select.raw.root",
            "top1_select.raw.root", "top2_select.raw.root", "top3_select.raw.root"]
    indir = "/eos/user/y/yofeng/LowPU/Selection/Wmunu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU/NTupleMod/Wmunu/ntuples_0_1/"
    logsuffix = "muonNtupleMod_wm_13"
    for fname in fnames:
        GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix)

    ## W -> enu
    macro = "eleNtupleMod.C"
    fnames = ["data_select.root",
            "we0_select.root", "we1_select.root", "we2_select.root",
            "wx0_select.root", "wx1_select.root", "wx2_select.root",
            "zxx_select.root",
            "ww_select.root", "wz_select.root", "zz_select.root",
            "top1_select.root", "top2_select.root", "top3_select.root"]
    indir = "/eos/user/y/yofeng/LowPU/Selection/Wenu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU/NTupleMod/Wenu/ntuples_0_1/"
    logsuffix = "eleNtupleMod_we_13"
    for fname in fnames:
        GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix)
