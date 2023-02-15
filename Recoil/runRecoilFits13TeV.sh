#! /bin/bash

IN13=/eos/user/y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Zmumu_pT20/merged/
IN13W=/eos/user/y/yofeng/Ntuples_LowPU/13TeV/20230213/Selections/Wmunu/merged/
OUT=/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_9_4_19/src/MitEwk13TeV/Recoil

LUMI13=200.87
#need to update the lumi^
#LUMI=1 #fitting for the fractions

#--------------------------------------------------------------------------------------------------
#-            CENTRAL VALUES
# -------------------------------------------------------------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWm.C+\(\"${IN13W}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_2G\",${LUMI13},0,0\)
root -l -b -q fitRecoilWm.C+\(\"${IN13W}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# 13 TeV Z MC
root -l -b -q fitRecoilZmm.C+\(\"${IN13}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_2G\",${LUMI13},0,0\)
# 13 TeV Z Data
root -l -b -q fitRecoilZmm.C+\(\"${IN13}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_2G_bkg_fixRoch\",${LUMI13},0,0\)

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ETA-BINNED
# -------------------------------------------------------------------------------------------------
for i in {1..3}; do
  echo ${i}
  # 13 TeV W MC
  root -l -b -q fitRecoilWm.C+\(\"${IN13W}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
  root -l -b -q fitRecoilWm.C+\(\"${IN13W}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
  # 13 TeV Z MC
  root -l -b -q fitRecoilZmm.C+\(\"${IN13}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
  # 13 TeV Z Data
  root -l -b -q fitRecoilZmm.C+\(\"${IN13}\",\"data.root\",2,2,0,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Eta${i}\",${LUMI13},${i},0\)
done

###################################################################################################
#--------------------------------------------------------------------------------------------------
#-            ROOKEYS PDF
# -------------------------------------------------------------------------------------------------
# 13 TeV W MC
root -l -b -q fitRecoilWm.C+\(\"${IN13W}\",2,2,1,1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmpMC_PF_13TeV_Keys\",${LUMI13},0,1\)
root -l -b -q fitRecoilWm.C+\(\"${IN13W}\",2,2,1,-1,0,\"met\",\"metPhi\",\"met\",\"${OUT}/WmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
# # 13 TeV Z MC
root -l -b -q fitRecoilZmm.C+\(\"${IN13}\",\"zmm.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmMC_PF_13TeV_Keys\",${LUMI13},0,1\)
# 13 TeV Z Data
root -l -b -q fitRecoilZmm.C+\(\"${IN13}\",\"data.root\",2,2,1,\"met\",\"metPhi\",\"met\",\"${OUT}/ZmmData_PF_13TeV_Keys\",${LUMI13},0,1\)

###################################################################################################
###################################################################################################

#rm *.so *.d *.pcm
