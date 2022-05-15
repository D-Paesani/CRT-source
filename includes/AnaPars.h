#ifndef AnaPars_H
#define AnaPars_H

    const int isRun182 = 0;

    const int centerMode = 0;
    const int strontiumMode = 0;

    const Long64_t maxEvToProcess3p1 = 1e4;
    const Long64_t maxEvToProcess3p2 = 1e4;
    const Long64_t maxEvToProcess3 = 1e8;

    const TString lutPrefix3p = "../../data/calibration/luts_s3p/";
    const TString lutPrefix3 = "../../data/calibration/luts_s3/";
    const TString lutChEqName = "_chargEq";
    const TString lutChEqErrName = "_chargEqErr";
    const TString lutBarLenName = "_barLen";
    const TString lutTimeOffsName = "_timeOff";
    const TString lutPeakTimeOffsName = "_peakTimeOff";
    const TString lutZetaOffName = "_zetaOff";
    const TString lutTimeDiffName = "_timeDiff";
    const TString lutTimeDiffErrName = "_timeDiffErr";

    const TString cutGPrefix = "../../data/calibration/cutg/";
    const TString cutGFormat = "cut%i%i";
    const int enableCutG = 0;

    const int enableOfflineEq = 0;

    const int sideNum = 2;
    const int scintNum = 8;

    const float scintVp = 14.3;
    const float scintL = 160.0;
    const float scintW = 2.5;


    const double minQCut = 50;

    const double maxQCut = 3000;
    const double maxVpeak = 1800;
    const double maxChi2 = 200000;
    const double maxQSharing = 40;
    const double chEqReference = 450;

    const float qFrom = 0;
    const float qTo = 1000; //deve essere 1000
    const int qBins = 1000; // deve essere 2000


    

    
#endif
