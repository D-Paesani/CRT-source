#include "AnaPars.h"
#include "TString.h"

using namespace std;



int GetSide(int n) { return (int)((n+1)>scintNum); }
int GetScint(int n) { return (int)(n - (GetSide(n)==1)*scintNum); }
int GetChan(int iSd, int iSc) { return (int)(iSc + iSd*scintNum); }

void NamerMatrix(int hN, TString& hTag, TString& hTitleTag) {
    hTag = Form("_%d_%d", GetSide(hN), GetScint(hN)); 
    hTitleTag = Form("[sd%d][sc%d] ", GetSide(hN), GetScint(hN));
}

void NamerArray(int n, TString& hTag, TString& hTitleTag) { 
  hTag = Form("_%d", n); 
  hTitleTag = Form("[ch%d] ", n); 
};

void SkipProc(TH1* histObj, int histN, int& histSkipFlag) { histSkipFlag=1; };


