{
    if (gSystem->Getenv("CMSSW_VERSION")) {
        gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
        gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaUtils.so");
        gROOT->Macro("$CMSSW_BASE/src/MitEwk13TeV/Utils/CEffUser2D.cc++");
        gROOT->Macro("$CMSSW_BASE/src/MitEwk13TeV/Utils/PdfDiagonalizer.cc++");
    }
    // Show which process needs debugging
    gInterpreter->ProcessLine(".! ps |grep root.exe");
}
