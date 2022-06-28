"""
script to generate and submit condor jobs
"""
import os
import time
from collections import OrderedDict

def GenerateExecutable(macro, indir, outdir, fname, logsuffix, nsec = 1, ith = 0, is5TeV = False):
    logsuffix += "_" + str(fname) + "_" + str(nsec) + "_" + str(ith)

    pwd = os.environ['PWD'];
    logname = "log_13TeV_new" if not is5TeV else "log_5TeV_new"
    jobname = pwd + "/" + logname + "/run_" + logsuffix
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

        cmd = "root -l -q -b "
        cmd += macro
        cmd += '++O\(\\"' + outdir + '\\",\\"' + indir
        if not is5TeV:
            cmd += '\\",\\"13TeV\\",\\"' + fname + '.root\\"'
        else:
            cmd += '\\",\\"5TeV\\",\\"' + fname + '.root\\"'

        if nsec != 1:
            cmd += ',' + str(nsec) + ',' + str(ith) 

        cmd += '\)'
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
+JobFlavour = "longlunch"
queue 1\n
""".format(jobname = jobname, here = pwd)

    with open(jobname + ".condor", 'w') as outfile:
        outfile.write(job_desc)
        outfile.close()

    os.system("condor_submit " + jobname + ".condor")
    time.sleep(0.05)


if __name__ == "__main__":
    # W -> munu
    macro = "muonNtupleMod.C"
    fnames = ["data_select.root", 
            "wm0_select.raw.root", "wm1_select.raw.root", "wm2_select.raw.root", 
            "wx0_select.raw.root", "wx1_select.raw.root", "wx2_select.raw.root", 
            "zxx_select.raw.root",
            "ww_select.raw.root", "wz_select.raw.root", "zz_select.raw.root",
            "top1_select.raw.root", "top2_select.raw.root", "top3_select.raw.root"]
    #fnames = ["wm0_select.raw.root"]

    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['wm0_select.raw.root'] = 200
    njobs['wm1_select.raw.root'] = 50
    njobs['wm2_select.raw.root'] = 10
    njobs['wx0_select.raw.root'] = 20
    njobs['wx1_select.raw.root'] = 10
    njobs['zxx_select.raw.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU/Selection/Wmunu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU/NTupleModTest4/Wmunu/ntuples_0_1/"
    logsuffix = "muonNtupleMod_wm_13"
    if 0:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob)

    ## W -> enu
    macro = "eleNtupleMod.C"
    fnames = ["data_select.root",
            "we0_select.root", "we1_select.root", "we2_select.root",
            "wx0_select.root", "wx1_select.root", "wx2_select.root",
            "zxx_select.root",
            "ww_select.root", "wz_select.root", "zz_select.root",
            "top1_select.root", "top2_select.root", "top3_select.root"]
    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['we0_select.root'] = 100
    njobs['we1_select.root'] = 50
    njobs['we2_select.root'] = 10
    njobs['wx0_select.root'] = 20
    njobs['wx1_select.root'] = 10
    njobs['zxx_select.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU/Selection/Wenu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU/NTupleModTest4/Wenu/ntuples_0_1/"
    logsuffix = "eleNtupleMod_we_13"
    if 0:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob)

    macro = "ZmmNTupleMod.C"
    fnames = ["data_select.root",
            "zmm_select.raw.root",
            "wx0_select.raw.root", "wx1_select.raw.root", "wx2_select.raw.root",
            "zxx_select.raw.root",
            "ww_select.raw.root","zz_select.raw.root", "wz_select.raw.root",
            "top1_select.raw.root", "top2_select.raw.root", "top3_select.raw.root",
            ]
    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['zmm_select.raw.root'] = 20
    njobs['zxx_select.raw.root'] = 10
    njobs['top3_select.raw.root'] = 10
    njobs['wx1_select.raw.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU/Selection_pT20/Zmumu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU/NTupleModTest4/Zmumu/ntuples_0_1/"
    logsuffix = "ZmmNTupleMod_zmumu_13"
    if 1:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob)

    macro = "ZeeNTupleMod.C"
    fnames = ["data_select.root",
            "zee_select.root",
            "wx0_select.root", "wx1_select.root", "wx2_select.root",
            "zxx_select.root",
            "ww_select.root","zz_select.root", "wz_select.root",
            "top1_select.root", "top2_select.root", "top3_select.root",
            ]
    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['zee_select.root'] = 20
    njobs['zxx_select.root'] = 10
    njobs['top3_select.root'] = 10
    njobs['wx1_select.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU/Selection_pT20/Zee/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU/NTupleModTest4/Zee/ntuples_0_1/"
    logsuffix = "ZeeNTupleMod_zee_13"
    if 0:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob)

    #
    # 5TeV results
    #
    macro = "muonNtupleMod.C"
    fnames = ["data_select.root", 
            "wm_select.raw.root",  
            "wx_select.raw.root",  
            "zxx_select.raw.root",
            "ww_select.raw.root", "wz_select.raw.root", "zz_select.raw.root",
            "top_select.raw.root"]

    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['wm_select.raw.root'] = 200
    njobs['wx_select.raw.root'] = 20
    njobs['zxx_select.raw.root'] = 20
    njobs['top_select.raw.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU_5TeV/Selection/Wmunu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/NTupleModTest4/Wmunu/ntuples_0_1/"
    logsuffix = "muonNtupleMod_wm_5"
    if 0:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob, True)

    ## W -> enu
    macro = "eleNtupleMod.C"
    fnames = ["data_select.root",
            "we_select.root", 
            "wx_select.root", 
            "zxx_select.root",
            "ww_select.root", "wz_select.root", "zz_select.root",
            "top_select.root", ]
    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['we_select.root'] = 100
    njobs['wx_select.root'] = 20
    njobs['zxx_select.root'] = 20
    njobs['top_select.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU_5TeV/Selection/Wenu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/NTupleModTest4/Wenu/ntuples_0_1/"
    logsuffix = "eleNtupleMod_we_5"
    if 0:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob, True)

    macro = "ZmmNTupleMod.C"
    fnames = ["data_select.root",
            "zmm_select.raw.root",
            "wx_select.raw.root", 
            "zxx_select.raw.root",
            "ww_select.raw.root","zz_select.raw.root", "wz_select.raw.root",
            "top_select.raw.root", 
            ]
    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['zmm_select.raw.root'] = 20
    njobs['zxx_select.raw.root'] = 10
    njobs['top_select.raw.root'] = 10
    njobs['wx_select.raw.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU_5TeV/Selection_pT20/Zmumu/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/NTupleModTest4/Zmumu/ntuples_0_1/"
    logsuffix = "ZmmNTupleMod_zmumu_5"
    if 1:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob, True)

    macro = "ZeeNTupleMod.C"
    fnames = ["data_select.root",
            "zee_select.root",
            "wx_select.root", 
            "zxx_select.root",
            "ww_select.root","zz_select.root", "wz_select.root",
            "top_select.root", 
            ]
    njobs = OrderedDict()
    for fname in fnames:
        njobs[fname] = 3
    njobs['data_select.root'] = 10
    njobs['zee_select.root'] = 20
    njobs['zxx_select.root'] = 10
    njobs['top_select.root'] = 10
    njobs['wx_select.root'] = 10
    indir = "/eos/user/y/yofeng/LowPU_5TeV/Selection_pT20/Zee/ntuples_0_1/"
    outdir = "/eos/user/y/yofeng/LowPU_5TeV/NTupleModTest4/Zee/ntuples_0_1/"
    logsuffix = "ZeeNTupleMod_zee_5"
    if 0:
        for fname in fnames:
            for ijob in range(njobs[fname]):
                GenerateExecutable(macro, indir, outdir, fname.replace(".root", ""), logsuffix, njobs[fname], ijob, True)
