"""
script to generate and submit condor jobs
"""
import os
import time

def GenerateExecutable(cmd, logsuffix, ijob):
    pwd = os.environ['PWD'];
    jobname = pwd + "/log/run_" + logsuffix + "_" + str(ijob)
    
    with open(jobname + ".sh", "w") as bashscript:
        bashscript.write("#!/bin/bash\n")
        bashscript.write("\n")
        bashscript.write("cd /afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Recoil\n")
        bashscript.write("eval $(scramv1 runtime -sh);\n")
        bashscript.write("source SetEnv.sh")
        bashscript.write("\n")

        # make the temp directory, copy the macro and libraries there
        bashscript.write("# create temp directory and copy macro and config there\n")
        tmpname = "/tmp/" + logsuffix + "/" + str(ijob) + "/"
        bashscript.write("mkdir -p " + tmpname + "\n")
        bashscript.write("cp rootlogon.C " + tmpname + "\n")
        bashscript.write("cp fitRecoilZll.C " + tmpname + "\n")
        bashscript.write("cp fitRecoilWl.C " + tmpname + "\n")
        bashscript.write("cd " + tmpname + "\n")
        bashscript.write("\n")

        bashscript.write(cmd)

    # change it to executable
    os.system("chmod +x " + jobname + ".sh")

    job_desc = """Universe = vanilla
Executable = {jobname}.sh
Log        = {jobname}.log
getenv      = True
environment = "LS_SUBCWD={here}"
+JobFlavour = "longlunch"
queue 1\n
""".format(jobname = jobname, here = pwd)

    with open(jobname + ".condor", 'w') as outfile:
        outfile.write(job_desc)
        outfile.close()

    os.system("condor_submit " + jobname + ".condor")
    time.sleep(0.1)
    
    
def GetNJobs(ifilename):
    with open(ifilename, "r") as ifile:
        nproc = 0
        for line in ifile:
            if line.startswith("root "):
                nproc += 1
    return nproc

def ReadBashScript(ifilename, ijob):
    with open(ifilename, "r") as ifile:
        iproc = 0
        for line in ifile:
            line = line.strip()
            if line.startswith("root "):
                if iproc == ijob:
                    return line
                iproc += 1
                
    raise RuntimeError("Cannot find job %d in %s" % (ijob, ifilename))
                

if __name__ == "__main__":
    #
    # 13 TeV
    #
    jobscript = "runRecoilFits13TeV.sh"    
    njobs = GetNJobs(jobscript)
    print("Found %d jobs in %s" % (njobs, jobscript))
    
    for ijob in range(njobs):
        cmd = ReadBashScript(jobscript, ijob)
        print(cmd)
        GenerateExecutable(cmd, "13TeV", ijob) 
        
    #
    # 5 TeV
    #
    jobscript = "runRecoilFits5TeV.sh"
    njobs = GetNJobs(jobscript)
    print("Found %d jobs in %s" % (njobs, jobscript))
    
    for ijob in range(njobs):
        cmd = ReadBashScript(jobscript, ijob)
        print(cmd)
        GenerateExecutable(cmd, "5TeV", ijob)