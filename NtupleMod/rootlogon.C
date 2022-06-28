{
    if (gSystem->Getenv("CMSSW_VERSION")) {
        gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libBaconAnaDataFormats.so");
        gSystem->Load("$CMSSW_BASE/src/MitEwk13TeV/Utils/CEffUser2D_cc.so");
        gSystem->Load("$CMSSW_BASE/src/MitEwk13TeV/Utils/PdfDiagonalizer_cc.so");
    }
    // Show which process needs debugging
    gInterpreter->ProcessLine(".! ps |grep root.exe");
}
