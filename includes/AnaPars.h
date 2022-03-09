#ifndef AnaPars_H
#define AnaPars_H


    const int centerMode = 0;

    const Long64_t maxEvToProcess3p1 = 100e3;
    const Long64_t maxEvToProcess3p2 = 100e3;
    const Long64_t maxEvToProcess3 = 1e9;

    const TString lutPrefix3p = "./data/calibration/luts_s3p/";
    const TString lutPrefix3 = "./data/calibration/luts_s3/";
    const TString lutChEqName = "_chargEq";
    const TString lutBarLenName = "_barLen";
    const TString lutTimeOffsName = "_timeOff";
    const TString lutZetaOffName = "_zetaOff";

    const TString cutGPrefix = "./data/calibration/cutg/";
    const TString cutGFormat = "cut%i%i";
    const int enableCutG = 0;

    const int enableOfflineEq = 0;







    
    const int sideNum = 2;
    const int scintNum = 8;

    const float scintVp = 12.5;
    const float scintL = 160.0;

    const double minQCut = 200;
    const double maxQCut = 100000;
    const double maxVpeak = 1800;
    const double maxChi2 = 200000;
    const double maxQSharing = 40;
    const double chEqReference = 450;

    const float qFrom = 50;
    const float qTo = 1250;
    const int qBins = 200;






    

    
#endif
