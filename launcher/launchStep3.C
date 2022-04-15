
#include "../CRT_step3.C"

#define inFile_f "../../data/step2/%s_s2.root"
#define outFile_f "../../data/step3/%s_s3.root"

void launchStep3(TString run_name, TString mod_select = "") {

  gErrorIgnoreLevel = kFatal; //kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal

  if (mod_select != "T" && mod_select != "B" && mod_select != "") {cout<<"modSel must be empty, 'T' or 'B' !!"<<endl; return;}
  //opzione "" per legacy

  TString inFileName = Form(inFile_f, run_name.Data());
  TString runName = mod_select == "" ? run_name : run_name + "_" + mod_select;
  TString calibName = runName;
  TString outFileName = Form(outFile_f, runName.Data());


  TFile *fileOut = new TFile(outFileName, "RECREATE");

  cout<<endl<<"------------> Launching step3:"<<endl;
  cout<<"----> runName : "<<runName<<endl;
  cout<<"----> Input file: "<<inFileName<<endl;
  cout<<"----> Output file: "<<outFileName<<endl;
  cout<<"----> Calibration files: "<<calibName<<endl;
  cout<<"----> Module selector: "<<mod_select<<endl<<endl;


  Analysis *a = new Analysis(inFileName, fileOut, runName, calibName, mod_select);

  a->Loop();



}










