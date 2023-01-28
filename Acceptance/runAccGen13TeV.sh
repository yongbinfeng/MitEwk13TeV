#! /bin/bash

OUTDIR=${PWD}/ResultsGEN
FRAC=0.5


########################################
##             W->munu
########################################

# # ################################################################################
# ##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -b -q computeAccGenWl.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmp_13TeV_amcPythia_dressed\",\"ACC_wmp_13_amcPythia_dressed\",1,1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmm_13TeV_amcPythia_dressed\",\"ACC_wmm_13_amcPythia_dressed\",1,-1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wm_13TeV_amcPythia_dressed\",\"ACC_wm_13_amcPythia_dressed\",1,0,0,1,${FRAC}\)

# # ################################################################################
# ##########     aMCnlo+Pythia with pT reweight for RESUMMATION
#root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmp_13TeV_amcPythia_ptWeight\",\"ACC_wmp_13_amcPythia_ptWeight\",1,1\)
#root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wmm_13TeV_amcPythia_ptWeight\",\"ACC_wmm_13_amcPythia_ptWeight\",1,-1\)
#root -l -q computeAccGenWm.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wm_13TeV_amcPythia_ptWeight\",\"ACC_wm_13_amcPythia_ptWeight\",1,0\)

# # ################################################################################
# ##########     Powheg+Pythia for RESUMMATION
root -l -b -q computeAccGenWl.C+\(\"wm13_powheg.conf\",\"${OUTDIR}/GEN_wmp_13TeV_powPythia_dressed\",\"ACC_wmp_13_powPythia_dressed\",1,1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"wm13_powheg.conf\",\"${OUTDIR}/GEN_wmm_13TeV_powPythia_dressed\",\"ACC_wmm_13_powPythia_dressed\",1,-1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"wm13_powheg.conf\",\"${OUTDIR}/GEN_wm_13TeV_powPythia_dressed\",\"ACC_wm_13_powPythia_dressed\",1,0,0,1,${FRAC}\)

# # ################################################################################
# ##########     Powheg+Photos for FSR
root -l -b -q computeAccGenWl.C+\(\"wm13_photos.conf\",\"${OUTDIR}/GEN_wmp_13TeV_powPhotos_dressed\",\"ACC_wmp_13_powPhotos_dressed\",1,1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"wm13_photos.conf\",\"${OUTDIR}/GEN_wmm_13TeV_powPhotos_dressed\",\"ACC_wmm_13_powPhotos_dressed\",1,-1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"wm13_photos.conf\",\"${OUTDIR}/GEN_wm_13TeV_powPhotos_dressed\",\"ACC_wm_13_powPhotos_dressed\",1,0,0,1,${FRAC}\)


# # ################################################################################
# # ##########     Powheg+Pythia for EWK
root -l -b -q computeAccGenWl.C+\(\"wm13_pythia.conf\",\"${OUTDIR}/GEN_wmp_13TeV_powPythia_dressed\",\"ACC_wmp_13_powPythia_dressed\",1,1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"wm13_pythia.conf\",\"${OUTDIR}/GEN_wmm_13TeV_powPythia_dressed\",\"ACC_wmm_13_powPythia_dressed\",1,-1,0,1,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"wm13_pythia.conf\",\"${OUTDIR}/GEN_wm_13TeV_powPythia_dressed\",\"ACC_wm_13_powPythia_dressed\",1,0,0,1,${FRAC}\)


########################################
##             W->enu
########################################

# ################################################################################
##########     aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -b -q computeAccGenWl.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wep_13TeV_amcPythia_dressed\",\"ACC_wep_13_amcPythia_dressed\",1,1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wem_13TeV_amcPythia_dressed\",\"ACC_wem_13_amcPythia_dressed\",1,-1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"w13.conf\",\"${OUTDIR}/GEN_we_13TeV_amcPythia_dressed\",\"ACC_we_13_amcPythia_dressed\",1,0,0,0,${FRAC}\)

# ################################################################################
##########     aMCnlo+Pythia for RESUMMATION with pT weight
#root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wep_13TeV_amcPythia_ptWeight\",\"ACC_wep_13_amcPythia_ptWeight\",1,1\)
#root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_wem_13TeV_amcPythia_ptWeight\",\"ACC_wem_13_amcPythia_ptWeight\",1,-1\)
#root -l -q computeAccGenWe.C+\(\"w13.conf\",\"${OUTDIR}/GEN_we_13TeV_amcPythia_ptWeight\",\"ACC_we_13_amcPythia_ptWeight\",1,0\)

# ################################################################################
##########     Powheg+Pythia for RESUMMATION
root -l -b -q computeAccGenWl.C+\(\"we13_powheg.conf\",\"${OUTDIR}/GEN_wep_13TeV_powPythia_dressed\",\"ACC_wep_13_powPythia_dressed\",1,1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"we13_powheg.conf\",\"${OUTDIR}/GEN_wem_13TeV_powPythia_dressed\",\"ACC_wem_13_powPythia_dressed\",1,-1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"we13_powheg.conf\",\"${OUTDIR}/GEN_we_13TeV_powPythia_dressed\",\"ACC_we_13_powPythia_dressed\",1,0,0,0,${FRAC}\)

# ################################################################################
##########     Powheg+Photos for EWK
root -l -b -q computeAccGenWl.C+\(\"we13_photos.conf\",\"${OUTDIR}/GEN_wep_13TeV_powLOPhotos_dressed\",\"ACC_wep_13_powLOPhotos_dressed\",1,1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"we13_photos.conf\",\"${OUTDIR}/GEN_wem_13TeV_powLOPhotos_dressed\",\"ACC_wem_13_powLOPhotos_dressed\",1,-1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"we13_photos.conf\",\"${OUTDIR}/GEN_we_13TeV_powLOPhotos_dressed\",\"ACC_we_13_powLOPhotos_dressed\",1,0,0,0,${FRAC}\)

# # ################################################################################
# ##########     powheg+Pythia for EWK
root -l -b -q computeAccGenWl.C+\(\"we13_pythia.conf\",\"${OUTDIR}/GEN_wep_13TeV_powLOPythia_dressed\",\"ACC_wep_13_powLOPythia_dressed\",1,1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"we13_pythia.conf\",\"${OUTDIR}/GEN_wem_13TeV_powLOPythia_dressed\",\"ACC_wem_13_powLOPythia_dressed\",1,-1,0,0,${FRAC}\)
root -l -b -q computeAccGenWl.C+\(\"we13_pythia.conf\",\"${OUTDIR}/GEN_we_13TeV_powLOPythia_dressed\",\"ACC_we_13_powLOPythia_dressed\",1,0,0,0,${FRAC}\)


########################################
##             Z->mm
########################################

# ################################################################################
# #########             aMCnlo+Pythia for RESUMMATION, QCD, and PDF
root -l -b -q computeAccGenZll.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zmm_13TeV_amcPythia_dressed\",\"ACC_zmm_13_amcPythia_dressed\",1,0,1,${FRAC}\)
# ################################################################################
# #########             aMCnlo+Pythia for RESUMMATION with pT weights
#root -l -q computeAccGenZmm.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zmm_13TeV_amcPythia_ptWeight\",\"ACC_zmm_13_amcPythia_ptWeight\",1\)
# # ################################################################################
# # #########             Powheg+minlo+Pythia for RESUMMATION
root -l -b -q computeAccGenZll.C+\(\"zmm13_minlo.conf\",\"${OUTDIR}/GEN_zmm_13TeV_powPythia_dressed\",\"ACC_zmm_13_powPythia_dressed\",1,0,1,${FRAC}\)
# # ################################################################################
# # ##########            Powheg+Photos for EWK
root -l -b -q computeAccGenZll.C+\(\"zmm13_photos.conf\",\"${OUTDIR}/GEN_zmm_13TeV_powLOPhotos_dressed\",\"ACC_zmm_13_powLOPhotos_dressed\",1,0,1,${FRAC}\)
# # ################################################################################
# # #########             powheg+Pythia for EWK
root -l -b -q computeAccGenZll.C+\(\"zmm13_pythia.conf\",\"${OUTDIR}/GEN_zmm_13TeV_powLOPythia_dressed\",\"ACC_zmm_13_powLOPythia_dressed\",1,0,1,${FRAC}\)


########################################
##             Z->ee
########################################
# ################################################################################
# # ##########             aMCnlo+Pythia for RESUMMATION
root -l -b -q computeAccGenZll.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zee_13TeV_amcPythia_dressed\",\"ACC_zee_13_amcPythia_dressed\",1,0,0,${FRAC}\)
# ################################################################################
# # ##########             aMCnlo+Pythia for RESUMMATION
#root -l -q computeAccGenZee.C+\(\"z13.conf\",\"${OUTDIR}/GEN_zee_13TeV_amcPythia_ptWeight\",\"ACC_zee_13_amcPythia_ptWeight\",1\)
# ################################################################################
# # ##########             Powheg+minlo+Pythia for RESUMMATION
root -l -b -q computeAccGenZll.C+\(\"zee13_minlo.conf\",\"${OUTDIR}/GEN_zee_13TeV_powPythia_dressed\",\"ACC_zee_13_powPythia_dressed\",1,0,0,${FRAC}\)
################################################################################
# # ##########             Powheg+Photos for EWK
root -l -b -q computeAccGenZll.C+\(\"zee13_photos.conf\",\"${OUTDIR}/GEN_zee_13TeV_powLOPhotos_dressed\",\"ACC_zee_13_powLOPhotos_dressed\",1,0,0,${FRAC}\)
################################################################################
# ##########             powheg+Pythia for EWK
root -l -b -q computeAccGenZll.C+\(\"zee13_pythia.conf\",\"${OUTDIR}/GEN_zee_13TeV_powLOPythia_dressed\",\"ACC_zee_13_powLOPythia_dressed\",1,0,0,${FRAC}\)

# rm *.so *.d
