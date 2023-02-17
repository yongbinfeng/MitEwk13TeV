#! /bin/bash
source SetEnv.sh

#--------------------------------------------------------------------------------------------------
#-            CENTRAL VALUES
# -------------------------------------------------------------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_2G\",${LUMI13},0,0\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# 13 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# 13 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_2G\",${LUMI13},0,0\)

#-----------------------------------------------
# electron
#-------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_2G\",${LUMI13},0,0,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_2G\",${LUMI13},0,0,0,1\)
# 13 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_2G\",${LUMI13},0,0,0,1\)
# 13 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_2G\",${LUMI13},0,0,0,1\)

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ETA-BINNED
# -------------------------------------------------------------------------------------------------
# central region
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Eta1\",${LUMI13},1,0\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Eta1\",${LUMI13},1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Eta1\",${LUMI13},1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Eta1\",${LUMI13},1,0\)
# mid region
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Eta2\",${LUMI13},2,0\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Eta2\",${LUMI13},2,0\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Eta2\",${LUMI13},2,0\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Eta2\",${LUMI13},2,0\)
# forward region
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Eta3\",${LUMI13},3,0\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Eta3\",${LUMI13},3,0\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Eta3\",${LUMI13},3,0\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Eta3\",${LUMI13},3,0\)

#----------------------------------------
# electron
#----------------------------------------
# central region
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_Eta1\",${LUMI13},1,0,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_Eta1\",${LUMI13},1,0,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_Eta1\",${LUMI13},1,0,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_Eta1\",${LUMI13},1,0,0,1\)
# mid region
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_Eta2\",${LUMI13},2,0,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_Eta2\",${LUMI13},2,0,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_Eta2\",${LUMI13},2,0,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_Eta2\",${LUMI13},2,0,0,1\)
# forward region
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_Eta3\",${LUMI13},3,0,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_Eta3\",${LUMI13},3,0,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_Eta3\",${LUMI13},3,0,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_Eta3\",${LUMI13},3,0,0,1\)

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ROOKEYS PDF
# -------------------------------------------------------------------------------------------------
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Keys\",${LUMI13},0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Keys\",${LUMI13},0,1\)

#----------------------------------------
# electron
#----------------------------------------
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)

###################################################################################################

#rm *.so *.d *.pcm
