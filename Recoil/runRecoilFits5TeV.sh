#! /bin/bash
source SetEnv.sh

#--------------------------------------------------------------------------------------------------
#-            CENTRAL VALUES
# -------------------------------------------------------------------------------------------------
# 5 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_2G\",${LUMI5},0,0,1,0\)
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_2G\",${LUMI5},0,0,1,0\)
# 5 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_2G\",${LUMI5},0,0,1,0\)
# 5 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_2G\",${LUMI5},0,0,1,0\)

#-----------------------------------------------
# electron
#-------------------------------------------
# 5 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_5TeV_2G\",${LUMI5},0,0,1,1\)
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_5TeV_2G\",${LUMI5},0,0,1,1\)
# 5 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_5TeV_2G\",${LUMI5},0,0,1,1\)
# 5 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_5TeV_2G\",${LUMI5},0,0,1,1\)

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ETA-BINNED
# -------------------------------------------------------------------------------------------------
# central region
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_Eta1\",${LUMI5},1,0,1,0\)
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_Eta1\",${LUMI5},1,0,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_Eta1\",${LUMI5},1,0,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_Eta1\",${LUMI5},1,0,1,0\)
# mid region
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_Eta2\",${LUMI5},2,0,1,0\)
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_Eta2\",${LUMI5},2,0,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_Eta2\",${LUMI5},2,0,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_Eta2\",${LUMI5},2,0,1,0\)
# forward region
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_Eta3\",${LUMI5},3,0,1,0\)
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_Eta3\",${LUMI5},3,0,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_Eta3\",${LUMI5},3,0,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_Eta3\",${LUMI5},3,0,1,0\)

#----------------------------------------
# electron
#----------------------------------------
# central region
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_5TeV_Eta1\",${LUMI5},1,0,1,1\)
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_5TeV_Eta1\",${LUMI5},1,0,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_5TeV_Eta1\",${LUMI5},1,0,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_5TeV_Eta1\",${LUMI5},1,0,1,1\)
# mid region
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_5TeV_Eta2\",${LUMI5},2,0,1,1\)
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_5TeV_Eta2\",${LUMI5},2,0,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_5TeV_Eta2\",${LUMI5},2,0,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_5TeV_Eta2\",${LUMI5},2,0,1,1\)
# forward region
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_5TeV_Eta3\",${LUMI5},3,0,1,1\)
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_5TeV_Eta3\",${LUMI5},3,0,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_5TeV_Eta3\",${LUMI5},3,0,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_5TeV_Eta3\",${LUMI5},3,0,1,1\)

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ROOKEYS PDF
# -------------------------------------------------------------------------------------------------
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_5TeV_Keys\",${LUMI5},0,1,1,0\)
root -l -b -q fitRecoilWl.C+\(\"${In5_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_5TeV_Keys\",${LUMI5},0,1,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_5TeV_Keys\",${LUMI5},0,1,1,0\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_5TeV_Keys\",${LUMI5},0,1,1,0\)

#----------------------------------------
# electron
#----------------------------------------
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_5TeV_Keys\",${LUMI5},0,1,1,1\)
root -l -b -q fitRecoilWl.C+\(\"${In5_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_5TeV_Keys\",${LUMI5},0,1,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_5TeV_Keys\",${LUMI5},0,1,1,1\)
root -l -b -q fitRecoilZll.C+\(\"${In5_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_5TeV_Keys\",${LUMI5},0,1,1,1\)

###################################################################################################
#rm *.so *.d *.pcm
