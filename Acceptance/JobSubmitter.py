import os
import time

def GetCommands(cmdfile):
    """Reads a file of commands and returns a list of commands."""
    # environment variables
    commands_com = []
    # running comands
    commands_run = []
    with open(cmdfile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line == '':
                continue
            if line.startswith('root'):
                commands_run.append(line)
            else:
                commands_com.append(line)
    return commands_com, commands_run

def GenerateExecutable(cmd_run, cmd_env, outdir, idx):
    """
    Generate the executable file for condor and submission file
    """
    pwd = os.environ['PWD']
    jobname = pwd + "/log/run_" + str(idx)
    with open(jobname + ".sh", "w") as bashscript:
        bashscript.write("#!/bin/bash\n")
        bashscript.write('\n')
        bashscript.write('echo "Starting job on " `date` #Date/time of start of job\n')
        bashscript.write('echo "Running on: `uname -a`" #Condor job is running on this node\n')
        bashscript.write('echo "System software: `cat /etc/redhat-release`" #Operating System on that node\n')
        bashscript.write('mkdir tmp' + idx + '\n')
        bashscript.write('cd tmp' + idx + '\n')
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
        bashscript.write("cd $CMSSW_BASE/src/MitEwk13TeV/Acceptance\n")
        bashscript.write("\n")

        for cmd in cmd_env:
            bashscript.write(cmd+"\n")
        bashscript.write(cmd_run+"\n")

        bashscript.write("xrdcp -r ./ResultsGEN/* root://cmseos.fnal.gov/" + outdir + "/\n")

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
    time.sleep(0.01)

    return

if __name__ == '__main__':
    # get the commands
    commands_com, commands_run = GetCommands('runAccGen13TeV.sh')
    print("commands_com: ", commands_com)
    print("commands_run: ", commands_run)

    print(len(commands_run), " jobs to be submitted. ")
    print(len(commands_com), " jobs to be submitted. ")

    outdir = "/store/user/yofeng/LowPUResults/TestAccept"
    # creat the outdir first
    os.system('/usr/bin/eos root://cmseos.fnal.gov rm -r ' + outdir)
    os.system('/usr/bin/eos root://cmseos.fnal.gov mkdir -p ' + outdir)

    # generate the executable
    for ijob in xrange(len(commands_run)):
        GenerateExecutable(commands_run[ijob], commands_com, outdir, str(ijob))