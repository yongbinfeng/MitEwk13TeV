#! /bin/bash

In13_Zmm=/eos/user/y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Zmumu_pT20/merged/
In13_Wm=/eos/user/y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Wmunu/merged/
In13_Zee=/eos/user/y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Zee_pT20/merged/
In13_We=/eos/user/y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Wenu/merged/
OUT=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Recoil/Results

LUMI13=200.87
#need to update the lumi^
#LUMI=1 #fitting for the fractions

#--------------------------------------------------------------------------------------------------
#-            CENTRAL VALUES
# -------------------------------------------------------------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_2G\",${LUMI13},0,0\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# 13 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# 13 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_2G_bkg_fixRoch\",${LUMI13},0,0\)

#-----------------------------------------------
# electron
#-------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_2G\",${LUMI13},0,0,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_2G\",${LUMI13},0,0,0,1\)
# 13 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_2G\",${LUMI13},0,0,0,1\)
# 13 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_2G_bkg_fixRoch\",${LUMI13},0,0,0,1\)

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ETA-BINNED
# -------------------------------------------------------------------------------------------------
for i in {1..3}; do
  echo ${i}
  # 13 TeV W MC
  root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
  root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
  # 13 TeV Z MC
  root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
  # 13 TeV Z Data
  root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
done

#----------------------------------------
# electron
#----------------------------------------
for i in {1..3}; do
  echo ${i}
  # 13 TeV W MC
  root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0,0,1\)
  root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0,0,1\)
  # 13 TeV Z MC
  root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0,0,1\)
  # 13 TeV Z Data
  root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_Eta${i}\",${LUMI13},${i},0,0,1\)
done

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ROOKEYS PDF
# -------------------------------------------------------------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Keys\",${LUMI13},0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_Wm}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
# # 13 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
# 13 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In13_Zmm}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Keys\",${LUMI13},0,1\)

#----------------------------------------
# electron
#----------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WepMC_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)
root -l -b -q fitRecoilWl.C+\(\"${In13_We}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WemMC_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)
# 13 TeV Z MC
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"zee.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeMC_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)
# 13 TeV Z Data
root -l -b -q fitRecoilZll.C+\(\"${In13_Zee}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZeeData_PF_13TeV_Keys\",${LUMI13},0,1,0,1\)

###################################################################################################

#rm *.so *.d *.pcm
