#ifndef MipSelection_h
#define MipSelection_h

#include <iostream>
#include <numeric>
#include <algorithm>

#include "AnaPars.h"

#include "TCutG.h"
using namespace std;



class MipSelection {
public:
  TCutG *cuts[2][8];

  void load_2dcuts(){
    for(int iSd=0; iSd<2; iSd++){
      for(int iSc=0; iSc<scintNum; iSc++){
        TFile f(Form("data/calibration/2dcuts/cut%i%i.root", iSd, iSc));
        cuts[iSd][iSc] = (TCutG*)f.Get("CUTG");
      }
    }
  }

  int isCoincidence(double Q[], int isc) { return (int)( Q[isc] > minQCut &&  Q[scintNum + isc] > minQCut ); }

  int isChargeGood(double Q[], int isc) {
    return (int)( Q[scintNum + isc] > minQCut && Q[scintNum + isc] < maxQCut && Q[isc] > minQCut && Q[isc] < maxQCut ); 
  }

  int isScintSaturated(double V[], int isc) { return (int)( V[isc] > maxVpeak || V[scintNum + isc] > maxVpeak ); }

  int isSaturated(double V) { return (int)( V > maxVpeak ); }

  int isEventSaturated(double V[])   {
    double max = *(max_element(V, V + 2*scintNum )); 
    return (max > maxVpeak);
  }

  int isShared(double Q[], int isc) { 
    if (isc > 0)            { if ( Q[isc-1] > maxQSharing || Q[scintNum+isc-1] > maxQSharing) { return 1; } }
    if (isc < scintNum - 1) { if ( Q[isc+1] > maxQSharing || Q[scintNum+isc+1] > maxQSharing) { return 1; } }
    return 0;
  }

  int isX2Good(double chi2[], int isc) {
    return (int)( chi2[scintNum + isc] > 0 && chi2[scintNum + isc] < maxChi2 && chi2[isc] > 0 && chi2[isc] < maxChi2 ); 
  }

  int hitsPrecheck(int ncry, int iside[]) { 
    if (ncry < 2) { return 0; }
    int sum = 0;
    accumulate(iside, iside + ncry, sum);
    if (sum == 0 || sum == ncry) { return 0;} 
    return 1;
  }

  int hitPrecheck(int isc, int iscint[], int ncry) { return 1 ; return (int)( count(iscint, iscint+ncry, isc) > 1 ); } /////////

  int isZetaGood(double Z) { return (int)( 0.9*Z < scintL/2 && 0.9*Z > -scintL/2 ); }

  int isTimeGood(double T) { return 1;}

  int ismuon(double Q[], double Z, int iSc){
    int side0ok = cuts[0][iSc]->IsInside(Z, Q[iSc]);
    int side1ok = cuts[1][iSc]->IsInside(Z, Q[iSc+scintNum]);
    return side0ok && side1ok;
  }

  int applyMipCuts() { return 1;} //whole mip selection

};


#endif



