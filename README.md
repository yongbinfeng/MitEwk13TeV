## Introduction

This folder contains the various modules needed to go from BACON ntuples to W/Z production cross section measurements. 

This code base has a twiki at: http://www.cmsaf.mit.edu/twiki/bin/view/CmsHep/WZBosonCrossSection and each module has its
own README.txt included. The general intention is that the TWIKI page should be read first for an explanation of the 
experimental/scientific background the overall analysis and each individual module. The README.txt's explain the technical
minutia of how to actually run the code as well as guidance on where modifications will be needed to use this code on new
data. 

## Modules

The Bacon ntuples are produced with the BaconProd module and read using the BaconAna module.

The order of modules is roughly as follows:

Selection:	     bacon ntuple	-> flat ntuple
		     (e scale cor)

EleScale:	     zee ntuple		-> electron scale/res corrections
		     raw sel ntuple

Efficiency:	     flat ntuple	-> lepton efficiency root files

Acceptance:	     flat ntuple	-> acceptance numbers
		     lep. eff. files

Recoil:		     flat ntuple	-> recoil correction root files

Signal Extraction:   flat ntuple	-> signal yield numbers
		     rec. corr. files

Summary Plots:	     signal yields	-> cross sections and plots
		     acc. numbers
		     syst. uncer. %'s

## Tools

The Tools folder contains a number of useful macros I needed along the way including ones for merging root files, merging and doing simple logic with json files, etc. I also used tools available at https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiCalc from the official CMS Luminosity calculation tools.

## Recipe

To make the analysis-level ntuples, first need to run the selections, then run the scripts in NtupleMod to apply lepton scale factors and recoil corrections, together with the uncertainties.

## ToDos
- Figure out how to run condor jobs
- Fix the muon-jet cleaning for ecal prefiring
- Fix the electron isolation calculation
- Fix the lepton energy scale / rochester correction uncertainty
- Implement the muon prefiring
