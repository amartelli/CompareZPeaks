#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "geometryUtils.h"
#include "ZPeakFitUtils.h"
#include "PUReweighting.h"
#include "GetScaleCorrection.h"
#include "GetExtraSmearing.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMath.h"



//-----------------
// global variables

std::map<std::string,TH1F*> h_nPV_MC;
std::map<std::string,TH1F*> h_nPV_DA;

TFile* outFile;
std::string extension;

std::string GetRunLabel(const int& runId);

TH1F* ratioHisto(TH1F* h_num, TH1F* h_den);

void DrawNVtx(const std::string& category, const std::string& label, const int& rebin,
              const std::string& outFileName);





int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>> drawNVtx::usage: " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  
  //----------------------
  // Parse the config file
  
  parseConfigFile(argv[1]);
  
  std::string inputFilesDA = gConfigParser -> readStringOption("Input::inputFilesDA");
  std::string inputFilesMC = gConfigParser -> readStringOption("Input::inputFilesMC");
  
  extension  = gConfigParser -> readStringOption("Options::extension");
  int maxEntries = gConfigParser -> readIntOption("Options::maxEntries");
  std::string dataLabel = gConfigParser -> readStringOption("Options::dataLabel");
  
  std::string outFilePath = gConfigParser -> readStringOption("Output::outFilePath");
  
  
  
  //------------------
  // Set style options
  
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  
  
  //----------
  // Get trees
  
  std::cout << std::endl;
  
  TChain* ntu_MC = new TChain("selected");
  FillChain(ntu_MC,inputFilesMC);
  std::cout << ">>>   MC: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  
  TChain* ntu_DA = new TChain("selected");
  FillChain(ntu_DA,inputFilesDA);
  std::cout << ">>> DATA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;
  
  if( ntu_MC->GetEntries() == 0 || ntu_DA->GetEntries() == 0 )
  {
    std::cout << ">>> drawNVtx::Error: at least one file is empty" << std::endl; 
    return -1;
  }
  
  
  
  //---------------------
  // Set branch addresses
  
  int runId;
  int nPV;
  int nPU;
  bool HLTfire;
  int* eleID = new int[2];
   
  ntu_DA -> SetBranchStatus("*",0);
  ntu_DA -> SetBranchStatus("runNumber",1); ntu_DA -> SetBranchAddress("runNumber",&runId);
  ntu_DA -> SetBranchStatus("nPV",      1); ntu_DA -> SetBranchAddress("nPV",      &nPV);
  ntu_DA -> SetBranchStatus("HLTfire",  1); ntu_DA -> SetBranchAddress("HLTfire",  &HLTfire);
  ntu_DA -> SetBranchStatus("eleID",    1); ntu_DA -> SetBranchAddress("eleID",    eleID);
  
  ntu_MC -> SetBranchStatus("runNumber",1); ntu_MC -> SetBranchAddress("runNumber",&runId);
  ntu_MC -> SetBranchStatus("nPV",      1); ntu_MC -> SetBranchAddress("nPV",      &nPV);  
  ntu_MC -> SetBranchStatus("nPU",      1); ntu_MC -> SetBranchAddress("nPU",      &nPU);  
  ntu_MC -> SetBranchStatus("HLTfire",  1); ntu_MC -> SetBranchAddress("HLTfire",  &HLTfire);
  ntu_MC -> SetBranchStatus("eleID",    1); ntu_MC -> SetBranchAddress("eleID",    eleID);
  
  
  // Define histograms
  std::vector<std::string> categories;
  
  categories.push_back("all");
  categories.push_back("all_HLTfire");
  categories.push_back("all_HLTfire_eleID");
  categories.push_back("all_noPUWeight");
  categories.push_back("Run2012AB");
  categories.push_back("Run2012AB_HLTfire");
  categories.push_back("Run2012AB_HLTfire_eleID");
  categories.push_back("Run2012AB_noPUWeight");
  categories.push_back("Run2012C");
  categories.push_back("Run2012C_HLTfire");
  categories.push_back("Run2012C_HLTfire_eleID");
  categories.push_back("Run2012C_noPUWeight");
  categories.push_back("Run2012D");
  categories.push_back("Run2012D_HLTfire");
  categories.push_back("Run2012D_HLTfire_eleID");
  categories.push_back("Run2012D_noPUWeight");
  
  for(unsigned int i = 0; i < categories.size(); ++i)
  {
    std::string category = categories.at(i);
    
    std::string histoName = "h_nPV_MC_"+category;
    h_nPV_MC[category] = new TH1F(histoName.c_str(),"",50,0.5,50.5);
    h_nPV_MC[category] -> Sumw2();
    
    histoName = "h_nPV_DA_"+category;
    h_nPV_DA[category] = new TH1F(histoName.c_str(),"",50,0.5,50.5);
    h_nPV_DA[category] -> Sumw2();
  }
  
  
  
  // define arrays for KS test
  std::map<std::string,std::vector<double> > masses;
  std::map<std::string,std::map<double,std::vector<double> > > masses_c;
  
  
  
  // pileup reweighting for MC
  std::map<std::string, std::map<float,float> > PUWeights;
  std::string PUDir(getenv("COMMONUTILS"));
  
  if( dataLabel == "Winter2013" )
  {
    PUWeights["Run2012AB"] = *(ComputePUWeights((PUDir+"/data/pileup/PUWeights_DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent_Run2012AB.root").c_str(),"h_PUweights",false));
    PUWeights["Run2012C"]  = *(ComputePUWeights((PUDir+"/data/pileup/PUWeights_DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent_Run2012C.root").c_str(), "h_PUweights",false));
    PUWeights["Run2012D"]  = *(ComputePUWeights((PUDir+"/data/pileup/PUWeights_DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent_Run2012D.root").c_str(), "h_PUweights",false));
  }
  
  
  
  // Loop over entries
  std::cout << std::endl;
  std::cout << ">>> Read data from MC sample" << std::endl;
  
  int nEntries_MC = ntu_MC -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_MC; ++ientry)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   MC entry " << ientry << " / " << nEntries_MC << "\r"<< std::flush;
    ntu_MC -> GetEntry(ientry);
    
    
    // variables
    std::string runLabel = GetRunLabel(runId);
    float ww = (PUWeights[runLabel])[int(nPU+0.5)];
    bool passEleID = false;
    if( ((eleID[0] & 6) == 6) && ((eleID[1] & 6) == 6) ) passEleID = true;
    
    
    // selections
    //if( isZ != 1 ) continue;
    
    
    // fill histograms
    h_nPV_MC["all"]             ->  Fill( nPV,ww );
    if( HLTfire ) h_nPV_MC["all_HLTfire"]     ->  Fill( nPV,ww );
    if( HLTfire && passEleID ) h_nPV_MC["all_HLTfire_eleID"]     ->  Fill( nPV,ww );
    h_nPV_MC["all_noPUWeight"]  ->  Fill( nPV,1. );
    
    if( runLabel == "Run2012AB" )
    {
      h_nPV_MC["Run2012AB"]             ->  Fill( nPV,ww );
      if( HLTfire ) h_nPV_MC["Run2012AB_HLTfire"]     ->  Fill( nPV,ww );
      if( HLTfire && passEleID ) h_nPV_MC["Run2012AB_HLTfire_eleID"]     ->  Fill( nPV,ww );
      h_nPV_MC["Run2012AB_noPUWeight"]  ->  Fill( nPV,1. );
    }
    
    if( runLabel == "Run2012C" )
    {
      h_nPV_MC["Run2012C"]             ->  Fill( nPV,ww );
      if( HLTfire ) h_nPV_MC["Run2012C_HLTfire"]     ->  Fill( nPV,ww );
      if( HLTfire && passEleID ) h_nPV_MC["Run2012C_HLTfire_eleID"]     ->  Fill( nPV,ww );
      h_nPV_MC["Run2012C_noPUWeight"]  ->  Fill( nPV,1. );
    }
    
    if( runLabel == "Run2012D" )
    {
      h_nPV_MC["Run2012D"]             ->  Fill( nPV,ww );
      if( HLTfire ) h_nPV_MC["Run2012D_HLTfire"]     ->  Fill( nPV,ww );
      if( HLTfire && passEleID ) h_nPV_MC["Run2012D_HLTfire_eleID"]     ->  Fill( nPV,ww );
      h_nPV_MC["Run2012D_noPUWeight"]  ->  Fill( nPV,1. );
    }
  }
  std::cout << std::endl;
  
  
  int nEntries_DA = ntu_DA -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_DA; ++ientry)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   DA entry " << ientry << " / " << nEntries_DA << "\r"<< std::flush;
    ntu_DA -> GetEntry(ientry);
    
    
    // variables
    std::string runLabel = GetRunLabel(runId);
    float ww = 1.;
    bool passEleID = false;
    if( ((eleID[0] & 6) == 6) && ((eleID[1] & 6) == 6) ) passEleID = true;
    
    
    // selections
    //if( isZ != 1 ) continue;
    
    
    // fill histograms
    h_nPV_DA["all"]             ->  Fill( nPV,ww );
    if( HLTfire ) h_nPV_DA["all_HLTfire"]     ->  Fill( nPV,ww );
    if( HLTfire && passEleID ) h_nPV_DA["all_HLTfire_eleID"]     ->  Fill( nPV,ww );
    h_nPV_DA["all_noPUWeight"]  ->  Fill( nPV,1. );
    
    if( runLabel == "Run2012AB" )
    {
      h_nPV_DA["Run2012AB"]             ->  Fill( nPV,ww );
      if( HLTfire ) h_nPV_DA["Run2012AB_HLTfire"]     ->  Fill( nPV,ww );
      if( HLTfire && passEleID ) h_nPV_DA["Run2012AB_HLTfire_eleID"]     ->  Fill( nPV,ww );
      h_nPV_DA["Run2012AB_noPUWeight"]  ->  Fill( nPV,1. );
    }
    
    if( runLabel == "Run2012C" )
    {
      h_nPV_DA["Run2012C"]             ->  Fill( nPV,ww );
      if( HLTfire ) h_nPV_DA["Run2012C_HLTfire"]     ->  Fill( nPV,ww );
      if( HLTfire && passEleID ) h_nPV_DA["Run2012C_HLTfire_eleID"]     ->  Fill( nPV,ww );
      h_nPV_DA["Run2012C_noPUWeight"]  ->  Fill( nPV,1. );
    }
    
    if( runLabel == "Run2012D" )
    {
      h_nPV_DA["Run2012D"]             ->  Fill( nPV,ww );
      if( HLTfire ) h_nPV_DA["Run2012D_HLTfire"]     ->  Fill( nPV,ww );
      if( HLTfire && passEleID ) h_nPV_DA["Run2012D_HLTfire_eleID"]     ->  Fill( nPV,ww );
      h_nPV_DA["Run2012D_noPUWeight"]  ->  Fill( nPV,1. );
    }
  }
  std::cout << std::endl;
  
  
  
  // Drawings
  std::string folderName = outFilePath + "/" + dataLabel + "/";
  gSystem -> mkdir(folderName.c_str());
  
  std::string outFileName = folderName + "/drawNVtx";
  outFile = new TFile((outFileName+".root").c_str(),"RECREATE");
  
  TCanvas* dummy = new TCanvas("dummy","",0,0,700,600);
  dummy -> Print((outFileName+"."+extension+"[").c_str(),extension.c_str());
  
  DrawNVtx("all",              "[Run2012ABCD]",                 1,outFileName);
  DrawNVtx("all_HLTfire",      "[Run2012ABCD - with HLT]",      1,outFileName);
  DrawNVtx("all_HLTfire_eleID","[Run2012ABCD - with HLT+eleID]",1,outFileName);
  DrawNVtx("all_noPUWeight",   "[Run2012ABCD - no PU reweigh.]",1,outFileName);
  
  DrawNVtx("Run2012AB",              "[Run2012AB]",                 1,outFileName);
  DrawNVtx("Run2012AB_HLTfire",      "[Run2012AB - with HLT]",      1,outFileName);
  DrawNVtx("Run2012AB_HLTfire_eleID","[Run2012AB - with HLT+eleID]",1,outFileName);
  DrawNVtx("Run2012AB_noPUWeight",   "[Run2012AB - no PU reweigh.]",1,outFileName);
  
  DrawNVtx("Run2012C",              "[Run2012C]",                 1,outFileName);
  DrawNVtx("Run2012C_HLTfire",      "[Run2012C - with HLT]",      1,outFileName);
  DrawNVtx("Run2012C_HLTfire_eleID","[Run2012C - with HLT+eleID]",1,outFileName);
  DrawNVtx("Run2012C_noPUWeight",   "[Run2012C - no PU reweigh.]",1,outFileName);
  
  DrawNVtx("Run2012D",              "[Run2012D]",                 1,outFileName);
  DrawNVtx("Run2012D_HLTfire",      "[Run2012D - with HLT]",      1,outFileName);
  DrawNVtx("Run2012D_HLTfire_eleID","[Run2012D - with HLT+eleID]",1,outFileName);
  DrawNVtx("Run2012D_noPUWeight",   "[Run2012D - no PU reweigh.]",1,outFileName);
    
  outFile -> Close();
  dummy -> Print((outFileName+"."+extension+"]").c_str(),extension.c_str());
}






void DrawNVtx(const std::string& category, const std::string& label, const int& rebin,
              const std::string& outFileName)
{
  TCanvas* c = new TCanvas(("c_nPV_"+category).c_str(),("nPV - "+category).c_str(),0,0,700,600);
  c -> cd();
  
  TPad* p1 = new TPad("p1","p1",0., 0.30, 1., 1.);
  TPad* p2 = new TPad("p2","p2",0., 0., 1., 0.30);
  p1 -> SetLeftMargin(0.16);
  p1 -> SetRightMargin(0.13);
  p2 -> SetLeftMargin(0.16);
  p2 -> SetRightMargin(0.13);
  p1 -> SetTopMargin(0.08);
  p1 -> SetBottomMargin(0.02);
  p2 -> SetTopMargin(0.04);
  p2 -> SetBottomMargin(0.4);
  p1 -> Draw();
  p2 -> Draw();
  
  p1 -> cd();
  
  p1 -> SetGridx();
  p2 -> SetGridy();
  
  
  
  
  //-----------
  // draw peaks
  std::cout << ">>> drawNVtx::DrawNVtx::draw nPV distribution for category " << category << std::endl;
  
  if( h_nPV_MC[category]->Integral() > 0 )
    h_nPV_MC[category] -> Scale( 1. * h_nPV_DA[category]->Integral() / h_nPV_MC[category]->Integral() );
  
  char axisTitle[50];
  h_nPV_MC[category] -> Rebin(rebin);
  sprintf(axisTitle,"events / %.2e",h_nPV_MC[category]->GetBinWidth(1));
  h_nPV_MC[category] -> GetXaxis() -> SetRangeUser(75.,104.999);
  //h_nPV_MC[category] -> GetXaxis() -> SetLabelSize(0.04);
  h_nPV_MC[category] -> GetXaxis() -> SetLabelSize(0.);
  h_nPV_MC[category] -> GetXaxis() -> SetLabelFont(42);
  //h_nPV_MC[category] -> GetXaxis() -> SetTitleSize(0.05);
  h_nPV_MC[category] -> GetXaxis() -> SetTitleSize(0.);
  h_nPV_MC[category] -> GetXaxis() -> SetTitleOffset(1.20);
  h_nPV_MC[category] -> GetXaxis() -> SetTitle(("N_{PV}   "+label).c_str());
  //h_nPV_MC[category] -> GetYaxis() -> SetLabelSize(0.04);
  h_nPV_MC[category] -> GetYaxis() -> SetLabelSize(0.057);
  h_nPV_MC[category] -> GetYaxis() -> SetLabelFont(42);
  //h_nPV_MC[category] -> GetYaxis() -> SetTitleSize(0.05);
  h_nPV_MC[category] -> GetYaxis() -> SetTitleSize(0.071);
  h_nPV_MC[category] -> GetYaxis() -> SetTitleOffset(1.22);
  h_nPV_MC[category] -> GetYaxis() -> SetTitle(axisTitle);
  
  h_nPV_MC[category] -> SetLineWidth(1);
  h_nPV_MC[category] -> SetLineColor(kRed);
  h_nPV_MC[category] -> SetFillColor(kYellow);
  h_nPV_MC[category] -> SetMarkerColor(kRed);
  h_nPV_MC[category] -> SetMarkerSize(0);
  gPad->Update();
  
  h_nPV_DA[category] -> Rebin(rebin);
  sprintf(axisTitle,"events / %.2e",h_nPV_DA[category]->GetBinWidth(1));
  h_nPV_DA[category] -> GetXaxis() -> SetRangeUser(75.,104.999);
  //h_nPV_DA[category] -> GetXaxis() -> SetLabelSize(0.04);
  h_nPV_DA[category] -> GetXaxis() -> SetLabelSize(0.);
  h_nPV_DA[category] -> GetXaxis() -> SetLabelFont(42);
  //h_nPV_DA[category] -> GetXaxis() -> SetTitleSize(0.05);
  h_nPV_DA[category] -> GetXaxis() -> SetTitleSize(0.);
  h_nPV_DA[category] -> GetXaxis() -> SetTitleOffset(1.20);  
  h_nPV_DA[category] -> GetXaxis() -> SetTitle(("N_{PV}   "+label).c_str());
  //h_nPV_DA[category] -> GetYaxis() -> SetLabelSize(0.04);
  h_nPV_DA[category] -> GetYaxis() -> SetLabelSize(0.057);
  h_nPV_DA[category] -> GetYaxis() -> SetLabelFont(42);
  //h_nPV_DA[category] -> GetYaxis() -> SetTitleSize(0.05);
  h_nPV_DA[category] -> GetYaxis() -> SetTitleSize(0.071);
  h_nPV_DA[category] -> GetYaxis() -> SetTitleOffset(1.22);  
  h_nPV_DA[category] -> GetYaxis() -> SetTitle(axisTitle);
  
  h_nPV_DA[category] -> SetLineWidth(1);
  h_nPV_DA[category] -> SetLineColor(kBlack);
  h_nPV_DA[category] -> SetMarkerColor(kBlack);
  h_nPV_DA[category] -> SetMarkerStyle(20);
  h_nPV_DA[category] -> SetMarkerSize(0.7);
  gPad->Update();
  
  
  // print the plots
  p1 -> cd();
  
  float maximum = h_nPV_MC[category] -> GetMaximum();
  if( h_nPV_DA[category]->GetMaximum() > maximum) maximum = h_nPV_DA[category]->GetMaximum();
  
  h_nPV_MC[category] -> SetMinimum(0.);
  h_nPV_MC[category] -> SetMaximum(1.05*maximum);
  h_nPV_MC[category] -> Draw("hist");
  h_nPV_MC[category] -> Draw("P,same");
  h_nPV_DA[category] -> Draw("P,sames");
  
  
  
  p2 -> cd();
  
  TH1F* ratio_MC = ratioHisto(h_nPV_MC[category],h_nPV_MC[category]);
  TH1F* ratio_DA = ratioHisto(h_nPV_DA[category],h_nPV_MC[category]);
  
  ratio_MC -> SetFillColor(ratio_MC->GetLineColor());
  ratio_MC -> SetFillStyle(3001);
  ratio_MC -> Draw("E3");
  ratio_DA -> Draw("same");
  
  gPad -> Update();
  
  
  
  c -> Print((outFileName+"."+extension).c_str(),extension.c_str());
  delete c;
  
  outFile -> cd();
  
  h_nPV_MC[category] -> Write();
  h_nPV_DA[category] -> Write();
  
}



TH1F* ratioHisto(TH1F* h_num, TH1F* h_den)
{
  TH1F* h_ratio = (TH1F*)( h_num->Clone() );
  h_ratio -> Divide(h_den);
  
  double yMin = 0.75;
  double yMax = 1.25;
  h_ratio -> GetYaxis() -> SetRangeUser(yMin,yMax);
  h_ratio -> GetXaxis() -> SetLabelSize(0.13);
  h_ratio -> GetYaxis() -> SetLabelSize(0.13);
  h_ratio -> GetXaxis() -> SetLabelFont(42);
  h_ratio -> GetYaxis() -> SetLabelFont(42);
  h_ratio -> GetXaxis() -> SetTitleSize(0.17);
  h_ratio -> GetYaxis() -> SetTitleSize(0.17);
  h_ratio -> GetYaxis() -> SetTitleFont(42);
  h_ratio -> GetXaxis() -> SetTitleOffset(1.00);
  h_ratio -> GetYaxis() -> SetTitleOffset(0.50);
  h_ratio -> GetYaxis() -> SetNdivisions(204);
  h_ratio -> GetXaxis() -> SetTitle(h_num->GetXaxis()->GetTitle());
  h_ratio -> GetYaxis() -> SetTitle("ratio");
  
  return h_ratio;
}



std::string GetRunLabel(const int& runId)
{
  std::string label = "";
  
  if( runId >= 190456 && runId <= 196531 ) label = "Run2012AB";
  if( runId >= 198022 && runId <= 203742 ) label = "Run2012C";
  if( runId >= 203777 && runId <= 208686 ) label = "Run2012D";
  
  return label;
}
