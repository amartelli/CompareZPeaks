#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "geometryUtils.h"
#include "ZPeakFitUtils.h"
#include "PUReweighting.h"
#include "GetScaleCorrections.h"
#include "GetSmearings.h"

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
// global functions

void DrawZPeak(const std::string& category, const std::string& label, const int& rebin,
               const bool& doVoigtianFit, const bool& doCrystalBallFit,
               const std::string& EBEE, const std::string& outFileName);

TH1F* ratioHisto(TH1F* h_num, TH1F* h_den, const double& xMax = -999., TArrow** line = NULL);



//-----------------
// global variables

std::map<std::string,TH1F*> h_mZ_MC;
std::map<std::string,TH1F*> h_mZ_DA;

TFile* outFile;
std::string extension;






int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>> compareZPeaks::usage: " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  
  //----------------------
  // Parse the config file
  
  parseConfigFile(argv[1]);
  
  std::string inputFilesDA = gConfigParser -> readStringOption("Input::inputFilesDA");
  std::string inputFilesMC = gConfigParser -> readStringOption("Input::inputFilesMC");
  
  extension  = gConfigParser -> readStringOption("Options::extension");
  int maxEntries = gConfigParser -> readIntOption("Options::maxEntries");
  std::string enCorrType = gConfigParser -> readStringOption("Options::enCorrType");
  std::string dataLabel  = gConfigParser -> readStringOption("Options::dataLabel");
  
  bool applyPUWeight        = gConfigParser -> readBoolOption("Options::applyPUWeight");
  bool applyEnergyScaleCorr = gConfigParser -> readBoolOption("Options::applyEnergyScaleCorr");
  bool applyEnergySmearing  = gConfigParser -> readBoolOption("Options::applyEnergySmearing");
  
  std::string energyScaleCorrType = gConfigParser -> readStringOption("Options::energyScaleCorrType");
  std::string energySmearingType  = gConfigParser -> readStringOption("Options::energySmearingType");
  
  bool doVoigtianFit    = gConfigParser -> readBoolOption("Options::doVoigtianFit");
  bool doCrystalBallFit = gConfigParser -> readBoolOption("Options::doCrystalBallFit");
  
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
  
  TChain* ntu_MC = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
  FillChain(ntu_MC,inputFilesMC);
  std::cout << ">>>   MC: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  
  TChain* ntu_DA = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
  FillChain(ntu_DA,inputFilesDA);
  std::cout << ">>> DATA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;
  
  if( ntu_MC->GetEntries() == 0 || ntu_DA->GetEntries() == 0 )
  {
    std::cout << ">>> compareZPeaks::Error: at least one file is empty" << std::endl; 
    return -1;
  }
  
  
  
  //---------------------
  // Set branch addresses
  
  int runId;
  int isZ;
  int timeStampHigh;
  float mZ;
  int isEB1,isEB2;
  float scEta1,scEta2;
  float scERaw1,scERaw2;
  float scE1,scE2;
  float scEReg1,scEReg2;
  float E3x31,E3x32;
  float R91,R92;
  int seedIeta1,seedIeta2;
  int seedIphi1,seedIphi2;
  int seedIx1,seedIx2;
  int seedIy1,seedIy2;
  int seedIz1,seedIz2;
  float scLaserCorr1,scLaserCorr2;
  float seedLaserAlpha1,seedLaserAlpha2;
  float npu;
  
  ntu_DA -> SetBranchStatus("*",0);
  ntu_DA -> SetBranchStatus("runId",              1);   ntu_DA -> SetBranchAddress("runId",&runId);
  ntu_DA -> SetBranchStatus("isZ",                1);   ntu_DA -> SetBranchAddress("isZ",&isZ);
  ntu_DA -> SetBranchStatus("timeStampHigh",      1);   ntu_DA -> SetBranchAddress("timeStampHigh",&timeStampHigh);
  ntu_DA -> SetBranchStatus("ele1ele2_scM",       1);   ntu_DA -> SetBranchAddress("ele1ele2_scM",&mZ);
  ntu_DA -> SetBranchStatus("ele1_isEB",          1);   ntu_DA -> SetBranchAddress("ele1_isEB",&isEB1);
  ntu_DA -> SetBranchStatus("ele2_isEB",          1);   ntu_DA -> SetBranchAddress("ele2_isEB",&isEB2);
  ntu_DA -> SetBranchStatus("ele1_scEta",         1);   ntu_DA -> SetBranchAddress("ele1_scEta",&scEta1);
  ntu_DA -> SetBranchStatus("ele2_scEta",         1);   ntu_DA -> SetBranchAddress("ele2_scEta",&scEta2);
  ntu_DA -> SetBranchStatus("ele1_scERaw",        1);   ntu_DA -> SetBranchAddress("ele1_scERaw",&scERaw1);
  ntu_DA -> SetBranchStatus("ele2_scERaw",        1);   ntu_DA -> SetBranchAddress("ele2_scERaw",&scERaw2);
  ntu_DA -> SetBranchStatus("ele1_scE",           1);   ntu_DA -> SetBranchAddress("ele1_scE",&scE1);
  ntu_DA -> SetBranchStatus("ele2_scE",           1);   ntu_DA -> SetBranchAddress("ele2_scE",&scE2);
  if( enCorrType == "stdSC" )
  {
    ntu_DA -> SetBranchStatus("ele1_scE",1);   ntu_DA -> SetBranchAddress("ele1_scE",&scEReg1);
    ntu_DA -> SetBranchStatus("ele2_scE",1);   ntu_DA -> SetBranchAddress("ele2_scE",&scEReg2);
  }
  if( enCorrType == "eleTunedReg" )
  {
    ntu_DA -> SetBranchStatus("ele1_scE_regression",1);   ntu_DA -> SetBranchAddress("ele1_scE_regression",&scEReg1);
    ntu_DA -> SetBranchStatus("ele2_scE_regression",1);   ntu_DA -> SetBranchAddress("ele2_scE_regression",&scEReg2);
  }
  if( enCorrType == "phoTunedReg" )
  {
    ntu_DA -> SetBranchStatus("ele1_scE_regression_PhotonTuned",1);   ntu_DA -> SetBranchAddress("ele1_scE_regression_PhotonTuned",&scEReg1);
    ntu_DA -> SetBranchStatus("ele2_scE_regression_PhotonTuned",1);   ntu_DA -> SetBranchAddress("ele2_scE_regression_PhotonTuned",&scEReg2);
  }
  ntu_DA -> SetBranchStatus("ele1_e3x3",          1);   ntu_DA -> SetBranchAddress("ele1_e3x3",&E3x31);
  ntu_DA -> SetBranchStatus("ele2_e3x3",          1);   ntu_DA -> SetBranchAddress("ele2_e3x3",&E3x32);
  ntu_DA -> SetBranchStatus("ele1_seedIeta",      1);   ntu_DA -> SetBranchAddress("ele1_seedIeta",&seedIeta1);
  ntu_DA -> SetBranchStatus("ele2_seedIeta",      1);   ntu_DA -> SetBranchAddress("ele2_seedIeta",&seedIeta2);
  ntu_DA -> SetBranchStatus("ele1_seedIphi",      1);   ntu_DA -> SetBranchAddress("ele1_seedIphi",&seedIphi1);
  ntu_DA -> SetBranchStatus("ele2_seedIphi",      1);   ntu_DA -> SetBranchAddress("ele2_seedIphi",&seedIphi2);
  ntu_DA -> SetBranchStatus("ele1_seedIx",        1);   ntu_DA -> SetBranchAddress("ele1_seedIx",&seedIx1);
  ntu_DA -> SetBranchStatus("ele2_seedIx",        1);   ntu_DA -> SetBranchAddress("ele2_seedIx",&seedIx2);
  ntu_DA -> SetBranchStatus("ele1_seedIy",        1);   ntu_DA -> SetBranchAddress("ele1_seedIy",&seedIy1);
  ntu_DA -> SetBranchStatus("ele2_seedIy",        1);   ntu_DA -> SetBranchAddress("ele2_seedIy",&seedIy2);
  ntu_DA -> SetBranchStatus("ele1_seedZside",     1);   ntu_DA -> SetBranchAddress("ele1_seedZside",&seedIz1);
  ntu_DA -> SetBranchStatus("ele2_seedZside",     1);   ntu_DA -> SetBranchAddress("ele2_seedZside",&seedIz2);
  ntu_DA -> SetBranchStatus("ele1_scLaserCorr",   1);   ntu_DA -> SetBranchAddress("ele1_scLaserCorr",&scLaserCorr1);
  ntu_DA -> SetBranchStatus("ele2_scLaserCorr",   1);   ntu_DA -> SetBranchAddress("ele2_scLaserCorr",&scLaserCorr2);
  ntu_DA -> SetBranchStatus("ele1_seedLaserAlpha",1);   ntu_DA -> SetBranchAddress("ele1_seedLaserAlpha",&seedLaserAlpha1);
  ntu_DA -> SetBranchStatus("ele2_seedLaserAlpha",1);   ntu_DA -> SetBranchAddress("ele2_seedLaserAlpha",&seedLaserAlpha2);
  
  ntu_MC -> SetBranchStatus("*",0);
  ntu_MC -> SetBranchStatus("isZ",                1);   ntu_MC -> SetBranchAddress("isZ",&isZ);
  ntu_MC -> SetBranchStatus("timeStampHigh",      1);   ntu_MC -> SetBranchAddress("timeStampHigh",&timeStampHigh);
  ntu_MC -> SetBranchStatus("ele1ele2_scM",       1);   ntu_MC -> SetBranchAddress("ele1ele2_scM",&mZ);
  ntu_MC -> SetBranchStatus("ele1_isEB",          1);   ntu_MC -> SetBranchAddress("ele1_isEB",&isEB1);
  ntu_MC -> SetBranchStatus("ele2_isEB",          1);   ntu_MC -> SetBranchAddress("ele2_isEB",&isEB2);
  ntu_MC -> SetBranchStatus("ele1_scEta",         1);   ntu_MC -> SetBranchAddress("ele1_scEta",&scEta1);
  ntu_MC -> SetBranchStatus("ele2_scEta",         1);   ntu_MC -> SetBranchAddress("ele2_scEta",&scEta2);
  ntu_MC -> SetBranchStatus("ele1_scERaw",        1);   ntu_MC -> SetBranchAddress("ele1_scERaw",&scERaw1);
  ntu_MC -> SetBranchStatus("ele2_scERaw",        1);   ntu_MC -> SetBranchAddress("ele2_scERaw",&scERaw2);
  ntu_MC -> SetBranchStatus("ele1_scE",           1);   ntu_MC -> SetBranchAddress("ele1_scE",&scE1);
  ntu_MC -> SetBranchStatus("ele2_scE",           1);   ntu_MC -> SetBranchAddress("ele2_scE",&scE2);
  if( enCorrType == "stdSC" )
  {
    ntu_MC -> SetBranchStatus("ele1_scE",1);   ntu_MC -> SetBranchAddress("ele1_scE",&scEReg1);
    ntu_MC -> SetBranchStatus("ele2_scE",1);   ntu_MC -> SetBranchAddress("ele2_scE",&scEReg2);
  }
  if( enCorrType == "eleTunedReg" )
  {
    ntu_MC -> SetBranchStatus("ele1_scE_regression",1);   ntu_MC -> SetBranchAddress("ele1_scE_regression",&scEReg1);
    ntu_MC -> SetBranchStatus("ele2_scE_regression",1);   ntu_MC -> SetBranchAddress("ele2_scE_regression",&scEReg2);
  }
  if( enCorrType == "phoTunedReg" )
  {
    ntu_MC -> SetBranchStatus("ele1_scE_regression_PhotonTuned",1);   ntu_MC -> SetBranchAddress("ele1_scE_regression_PhotonTuned",&scEReg1);
    ntu_MC -> SetBranchStatus("ele2_scE_regression_PhotonTuned",1);   ntu_MC -> SetBranchAddress("ele2_scE_regression_PhotonTuned",&scEReg2);
  }
  ntu_MC -> SetBranchStatus("ele1_e3x3",          1);   ntu_MC -> SetBranchAddress("ele1_e3x3",&E3x31);
  ntu_MC -> SetBranchStatus("ele2_e3x3",          1);   ntu_MC -> SetBranchAddress("ele2_e3x3",&E3x32);
  ntu_MC -> SetBranchStatus("ele1_seedIeta",      1);   ntu_MC -> SetBranchAddress("ele1_seedIeta",&seedIeta1);
  ntu_MC -> SetBranchStatus("ele2_seedIeta",      1);   ntu_MC -> SetBranchAddress("ele2_seedIeta",&seedIeta2);
  ntu_MC -> SetBranchStatus("ele1_seedIphi",      1);   ntu_MC -> SetBranchAddress("ele1_seedIphi",&seedIphi1);
  ntu_MC -> SetBranchStatus("ele2_seedIphi",      1);   ntu_MC -> SetBranchAddress("ele2_seedIphi",&seedIphi2);
  ntu_MC -> SetBranchStatus("ele1_seedIx",        1);   ntu_MC -> SetBranchAddress("ele1_seedIx",&seedIx1);
  ntu_MC -> SetBranchStatus("ele2_seedIx",        1);   ntu_MC -> SetBranchAddress("ele2_seedIx",&seedIx2);
  ntu_MC -> SetBranchStatus("ele1_seedIy",        1);   ntu_MC -> SetBranchAddress("ele1_seedIy",&seedIy1);
  ntu_MC -> SetBranchStatus("ele2_seedIy",        1);   ntu_MC -> SetBranchAddress("ele2_seedIy",&seedIy2);
  ntu_MC -> SetBranchStatus("ele1_seedZside",     1);   ntu_MC -> SetBranchAddress("ele1_seedZside",&seedIz1);
  ntu_MC -> SetBranchStatus("ele2_seedZside",     1);   ntu_MC -> SetBranchAddress("ele2_seedZside",&seedIz2);
  ntu_MC -> SetBranchStatus("ele1_scLaserCorr",   1);   ntu_MC -> SetBranchAddress("ele1_scLaserCorr",&scLaserCorr1);
  ntu_MC -> SetBranchStatus("ele2_scLaserCorr",   1);   ntu_MC -> SetBranchAddress("ele2_scLaserCorr",&scLaserCorr2);
  ntu_MC -> SetBranchStatus("ele1_seedLaserAlpha",1);   ntu_MC -> SetBranchAddress("ele1_seedLaserAlpha",&seedLaserAlpha1);
  ntu_MC -> SetBranchStatus("ele2_seedLaserAlpha",1);   ntu_MC -> SetBranchAddress("ele2_seedLaserAlpha",&seedLaserAlpha2);
  ntu_MC -> SetBranchStatus("PUit_TrueNumInteractions",1); ntu_MC -> SetBranchAddress("PUit_TrueNumInteractions",&npu);  
  
  
  // Define histograms
  std::vector<std::string> categories;
  
  categories.push_back("EB-EB");
  categories.push_back("EB-EB_hR9");
  categories.push_back("EB-EB_lR9");
  categories.push_back("EB-EB_eta>1");
  categories.push_back("EB-EB_eta<1");
  categories.push_back("EB-EB_hR9_eta<1");
  categories.push_back("EB-EB_lR9_eta<1");
  categories.push_back("EB-EB_hR9_eta>1");
  categories.push_back("EB-EB_lR9_eta>1");
  
  categories.push_back("EE-EE");
  categories.push_back("EE-EE_hR9");
  categories.push_back("EE-EE_lR9");
  categories.push_back("EE-EE_eta>2");
  categories.push_back("EE-EE_eta<2");
  categories.push_back("EE-EE_hR9_eta<2");
  categories.push_back("EE-EE_lR9_eta<2");
  categories.push_back("EE-EE_hR9_eta>2");
  categories.push_back("EE-EE_lR9_eta>2");
  
  for(unsigned int i = 0; i < categories.size(); ++i)
  {
    std::string category = categories.at(i);
    
    std::string histoName = "h_mZ_MC_"+category;
    h_mZ_MC[category] = new TH1F(histoName.c_str(),"",160,70.,110.);
    h_mZ_MC[category] -> Sumw2();
    
    histoName = "h_mZ_DA_"+category;
    h_mZ_DA[category] = new TH1F(histoName.c_str(),"",160,70.,110.);
    h_mZ_DA[category] -> Sumw2();
  }
  
  
  
  // define arrays for KS test
  std::map<std::string,std::vector<double> > masses;
  std::map<std::string,std::map<double,std::vector<double> > > masses_c;
  
  
  
  // pileup reweighting for MC
  std::map<float,float> PUWeights;
  std::string PUDir(getenv("COMMONUTILS"));
  
  if( dataLabel == "Moriond2013" )
    PUWeights = *(ComputePUweights(ntu_MC,(PUDir+"/pileup/pileup_69p3mb_true_Moriond2013.root").c_str(),false));
  if( dataLabel == "Winter2013" )
    PUWeights = *(ComputePUweights(ntu_MC,(PUDir+"/pileup/pileup_69p3mb_true_Winter2013.root").c_str(),false));
  
  
  
  // Loop over entries
  std::cout << std::endl;
  std::cout << ">>> Read data from MC sample" << std::endl;
  
  TRandom3 r;
  int nEntries_MC = ntu_MC -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_MC; ++ientry)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   MC entry " << ientry << " / " << nEntries_MC << "\r"<< std::flush;
    ntu_MC -> GetEntry(ientry);
    
    
    // variables
    R91 = E3x31/scERaw1;
    R92 = E3x32/scERaw2;
    
    if( applyEnergySmearing )
    {
      float energySmearing1 = gRandom->Gaus(1.,GetSmearings(scEta1,R91,dataLabel,energySmearingType));
      float energySmearing2 = gRandom->Gaus(1.,GetSmearings(scEta2,R92,dataLabel,energySmearingType));
      scEReg1 *= energySmearing1;
      scEReg2 *= energySmearing2;
    }
  
    float ww = 1.;
    if( applyPUWeight )
    {
      ww *= PUWeights[int(npu+0.5)];
    }  
    
    
    // selections
    if( isZ != 1 ) continue;
    if( fabs(scEta1) > 2.5    || fabs(scEta2) > 2.5  ) continue;
    if( fabs(scEta1) > 1.4442 && fabs(scEta1) < 1.56 ) continue;
    if( fabs(scEta2) > 1.4442 && fabs(scEta2) < 1.56 ) continue;
    
    
    // use regression energy
    mZ *= sqrt( scEReg1/scE1 * scEReg2/scE2 );
    
    
    // fill histograms
    if( (isEB1 == 1) && (isEB2 == 1) )
    {
      h_mZ_MC["EB-EB"]  ->  Fill( mZ,ww );
      
      if( (R91 > 0.94) &&
          (R92 > 0.94) )
      {
        h_mZ_MC["EB-EB_hR9"]  ->  Fill( mZ,ww );
      }
      
      if( (R91 < 0.94) &&
          (R92 < 0.94) )
      {
        h_mZ_MC["EB-EB_lR9"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) > 1.) &&
          (fabs(scEta2) > 1.) )
      {
        h_mZ_MC["EB-EB_eta>1"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) < 1.) &&
          (fabs(scEta2) < 1.) )
      {
        h_mZ_MC["EB-EB_eta<1"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) < 1. && R91 > 0.94) &&
          (fabs(scEta2) < 1. && R92 > 0.94) )
      {
        h_mZ_MC["EB-EB_hR9_eta<1"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) < 1. && R91 < 0.94) &&
          (fabs(scEta2) < 1. && R92 < 0.94) )
      {
        h_mZ_MC["EB-EB_lR9_eta<1"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) > 1. && R91 > 0.94) &&
          (fabs(scEta2) > 1. && R92 > 0.94) )
      {
        h_mZ_MC["EB-EB_hR9_eta>1"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) > 1. && R91 < 0.94) &&
          (fabs(scEta2) > 1. && R92 < 0.94) )
      {
        h_mZ_MC["EB-EB_lR9_eta>1"]  ->  Fill( mZ,ww );
      }
    }
    
    
    
    // fill histograms
    if( (isEB1 == 0) && (isEB2 == 0) )
    {
      h_mZ_MC["EE-EE"]  ->  Fill( mZ,ww );
      
      if( (R91 > 0.94) &&
          (R92 > 0.94) )
      {
        h_mZ_MC["EE-EE_hR9"]  ->  Fill( mZ,ww );
      }
      
      if( (R91 < 0.94) &&
          (R92 < 0.94) )
      {
        h_mZ_MC["EE-EE_lR9"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) > 2.) &&
          (fabs(scEta2) > 2.) )
      {
        h_mZ_MC["EE-EE_eta>2"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) < 2.) &&
          (fabs(scEta2) < 2.) )
      {
        h_mZ_MC["EE-EE_eta<2"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) < 2. && R91 > 0.94) &&
          (fabs(scEta2) < 2. && R92 > 0.94) )
      {
        h_mZ_MC["EE-EE_hR9_eta<2"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) < 2. && R91 < 0.94) &&
          (fabs(scEta2) < 2. && R92 < 0.94) )
      {
        h_mZ_MC["EE-EE_lR9_eta<2"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) > 2. && R91 > 0.94) &&
          (fabs(scEta2) > 2. && R92 > 0.94) )
      {
        h_mZ_MC["EE-EE_hR9_eta>2"]  ->  Fill( mZ,ww );
      }
      
      if( (fabs(scEta1) > 2. && R91 < 0.94) &&
          (fabs(scEta2) > 2. && R92 < 0.94) )
      {
        h_mZ_MC["EE-EE_lR9_eta>2"]  ->  Fill( mZ,ww );
      }
    }
  }
  std::cout << std::endl;
  
  
  
  std::cout << ">>> Read data from DATA sample" << std::endl;
  
  int nEntries_DA = ntu_DA -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_DA; ++ientry)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading DATA entry " << ientry << " / " << nEntries_DA << "\r" << std::flush;
    ntu_DA -> GetEntry(ientry);
    
    
    // variables
    R91 = E3x31/scERaw1;
    R92 = E3x32/scERaw2;
    
    if( applyEnergyScaleCorr )
    {
      scEReg1 *= GetScaleCorrections(scEta1,R91,runId,dataLabel,energyScaleCorrType);
      scEReg2 *= GetScaleCorrections(scEta2,R92,runId,dataLabel,energyScaleCorrType);
    }
    
    
    // selections
    if( isZ != 1 ) continue;
    if( fabs(scEta1) > 2.5000 || fabs(scEta2) > 2.5000 ) continue;
    if( fabs(scEta1) > 1.4442 && fabs(scEta1) < 1.5600 ) continue;
    if( fabs(scEta2) > 1.4442 && fabs(scEta2) < 1.5600 ) continue;
    
    
    // use regression energy
    mZ *= sqrt( scEReg1/scE1 * scEReg2/scE2 );
    
    
    // fill histograms
    if( (isEB1 == 1) && (isEB2 == 1) )
    {
      h_mZ_DA["EB-EB"]  ->  Fill( mZ );
      
      if( (R91 > 0.94) &&
          (R92 > 0.94) )
      {
        h_mZ_DA["EB-EB_hR9"]  ->  Fill( mZ );
      }
      
      if( (R91 < 0.94) &&
          (R92 < 0.94) )
      {
        h_mZ_DA["EB-EB_lR9"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 1.) &&
          (fabs(scEta2) > 1.) )
      {
        h_mZ_DA["EB-EB_eta>1"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 1.) &&
          (fabs(scEta2) < 1.) )
      {
        h_mZ_DA["EB-EB_eta<1"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 1. && R91 > 0.94) &&
          (fabs(scEta2) < 1. && R92 > 0.94) )
      {
        h_mZ_DA["EB-EB_hR9_eta<1"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 1. && R91 < 0.94) &&
          (fabs(scEta2) < 1. && R92 < 0.94) )
      {
        h_mZ_DA["EB-EB_lR9_eta<1"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 1. && R91 > 0.94) &&
          (fabs(scEta2) > 1. && R92 > 0.94) )
      {
        h_mZ_DA["EB-EB_hR9_eta>1"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 1. && R91 < 0.94) &&
          (fabs(scEta2) > 1. && R92 < 0.94) )
      {
        h_mZ_DA["EB-EB_lR9_eta>1"]  ->  Fill( mZ );
      }
    }
    
    
    
    // fill histograms
    if( (isEB1 == 0) && (isEB2 == 0) )
    {
      h_mZ_DA["EE-EE"]  ->  Fill( mZ );
      
      if( (R91 > 0.94) &&
          (R92 > 0.94) )
      {
        h_mZ_DA["EE-EE_hR9"]  ->  Fill( mZ );
      }
      
      if( (R91 < 0.94) &&
          (R92 < 0.94) )
      {
        h_mZ_DA["EE-EE_lR9"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 2.) &&
          (fabs(scEta2) > 2.) )
      {
        h_mZ_DA["EE-EE_eta>2"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 2.) &&
          (fabs(scEta2) < 2.) )
      {
        h_mZ_DA["EE-EE_eta<2"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 2. && R91 > 0.94) &&
          (fabs(scEta2) < 2. && R92 > 0.94) )
      {
        h_mZ_DA["EE-EE_hR9_eta<2"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 2. && R91 < 0.94) &&
          (fabs(scEta2) < 2. && R92 < 0.94) )
      {
        h_mZ_DA["EE-EE_lR9_eta<2"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 2. && R91 > 0.94) &&
          (fabs(scEta2) > 2. && R92 > 0.94) )
      {
        h_mZ_DA["EE-EE_hR9_eta>2"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 2. && R91 < 0.94) &&
          (fabs(scEta2) > 2. && R92 < 0.94) )
      {
        h_mZ_DA["EE-EE_lR9_eta>2"]  ->  Fill( mZ );
      }
    }
  }
  std::cout << std::endl;
  
  
  
  // Drawings
  std::string folderName = outFilePath + "/" + dataLabel + "/";
  gSystem -> mkdir(folderName.c_str());
  
  std::string outFileName = folderName + "/compareZPeaks";
  
  outFileName += "__" + enCorrType;
  if( applyPUWeight )        outFileName += "__MC_PUweight";
  if( applyEnergySmearing )  outFileName += "__MC_energySmearing_" + energySmearingType;
  if( applyEnergyScaleCorr ) outFileName += "__DA_energyScaleCorr_" + energyScaleCorrType;
  
  outFile = new TFile((outFileName+".root").c_str(),"RECREATE");
  
  TCanvas* dummy = new TCanvas("dummy","",0,0,700,600);
  dummy -> Print((outFileName+"."+extension+"[").c_str(),extension.c_str());
  
  DrawZPeak("EB-EB",          "[EB-EB]",                         1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_hR9",      "[EB-EB   R9 > 0.94]",             1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_lR9",      "[EB-EB   R9 < 0.94]",             1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_eta>1",    "[EB-EB   |#eta| > 1]",            1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_eta<1",    "[EB-EB   |#eta| < 1]",            1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_hR9_eta<1","[EB-EB   R9 > 0.94   |#eta| < 1]",1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_lR9_eta<1","[EB-EB   R9 < 0.94   |#eta| < 1]",1,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_hR9_eta>1","[EB-EB   R9 > 0.94   |#eta| > 1]",2,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  DrawZPeak("EB-EB_lR9_eta>1","[EB-EB   R9 < 0.94   |#eta| > 1]",2,doVoigtianFit,doCrystalBallFit,"EB",outFileName);
  
  DrawZPeak("EE-EE",          "[EE-EE]",                         1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_hR9",      "[EE-EE   R9 > 0.94]",             1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_lR9",      "[EE-EE   R9 < 0.94]",             1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_eta>2",    "[EE-EE   |#eta| > 2]",            1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_eta<2",    "[EE-EE   |#eta| < 2]",            1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_hR9_eta<2","[EE-EE   R9 > 0.94   |#eta| < 2]",2,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_lR9_eta<2","[EE-EE   R9 < 0.94   |#eta| < 2]",1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_hR9_eta>2","[EE-EE   R9 > 0.94   |#eta| > 2]",1,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  DrawZPeak("EE-EE_lR9_eta>2","[EE-EE   R9 < 0.94   |#eta| > 2]",2,doVoigtianFit,doCrystalBallFit,"EE",outFileName);
  
  outFile -> Close();
  dummy -> Print((outFileName+"."+extension+"]").c_str(),extension.c_str());
}






void DrawZPeak(const std::string& category, const std::string& label, const int& rebin,
               const bool& doVoigtianFit, const bool& doCrystalBallFit,
               const std::string& EBEE, const std::string& outFileName)
{
  TCanvas* c = new TCanvas(("c_mZ_"+category).c_str(),("mZ - "+category).c_str(),0,0,700,600);
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
  std::cout << ">>> compareZPeaks::DrawZPeak::draw peaks for category " << category << std::endl;
  
  if( h_mZ_MC[category]->Integral() > 0 )
    h_mZ_MC[category] -> Scale( 1. * h_mZ_DA[category]->Integral() / h_mZ_MC[category]->Integral() );
  
  char axisTitle[50];
  h_mZ_MC[category] -> Rebin(rebin);
  sprintf(axisTitle,"events / %.2e GeV/c^{2}",h_mZ_MC[category]->GetBinWidth(1));
  h_mZ_MC[category] -> GetXaxis() -> SetRangeUser(75.,104.999);
  //h_mZ_MC[category] -> GetXaxis() -> SetLabelSize(0.04);
  h_mZ_MC[category] -> GetXaxis() -> SetLabelSize(0.);
  h_mZ_MC[category] -> GetXaxis() -> SetLabelFont(42);
  //h_mZ_MC[category] -> GetXaxis() -> SetTitleSize(0.05);
  h_mZ_MC[category] -> GetXaxis() -> SetTitleSize(0.);
  h_mZ_MC[category] -> GetXaxis() -> SetTitleOffset(1.20);
  h_mZ_MC[category] -> GetXaxis() -> SetTitle(("m_{ee}   "+label).c_str());
  //h_mZ_MC[category] -> GetYaxis() -> SetLabelSize(0.04);
  h_mZ_MC[category] -> GetYaxis() -> SetLabelSize(0.057);
  h_mZ_MC[category] -> GetYaxis() -> SetLabelFont(42);
  //h_mZ_MC[category] -> GetYaxis() -> SetTitleSize(0.05);
  h_mZ_MC[category] -> GetYaxis() -> SetTitleSize(0.071);
  h_mZ_MC[category] -> GetYaxis() -> SetTitleOffset(1.22);
  h_mZ_MC[category] -> GetYaxis() -> SetTitle(axisTitle);
  
  h_mZ_MC[category] -> SetLineWidth(1);
  h_mZ_MC[category] -> SetLineColor(kRed);
  h_mZ_MC[category] -> SetFillColor(kYellow);
  h_mZ_MC[category] -> SetMarkerColor(kRed);
  h_mZ_MC[category] -> SetMarkerSize(0);
  gPad->Update();
  
  h_mZ_DA[category] -> Rebin(rebin);
  sprintf(axisTitle,"events / %.2e GeV/c^{2}",h_mZ_DA[category]->GetBinWidth(1));
  h_mZ_DA[category] -> GetXaxis() -> SetRangeUser(75.,104.999);
  //h_mZ_DA[category] -> GetXaxis() -> SetLabelSize(0.04);
  h_mZ_DA[category] -> GetXaxis() -> SetLabelSize(0.);
  h_mZ_DA[category] -> GetXaxis() -> SetLabelFont(42);
  //h_mZ_DA[category] -> GetXaxis() -> SetTitleSize(0.05);
  h_mZ_DA[category] -> GetXaxis() -> SetTitleSize(0.);
  h_mZ_DA[category] -> GetXaxis() -> SetTitleOffset(1.20);  
  h_mZ_DA[category] -> GetXaxis() -> SetTitle(("m_{ee}   "+label).c_str());
  //h_mZ_DA[category] -> GetYaxis() -> SetLabelSize(0.04);
  h_mZ_DA[category] -> GetYaxis() -> SetLabelSize(0.057);
  h_mZ_DA[category] -> GetYaxis() -> SetLabelFont(42);
  //h_mZ_DA[category] -> GetYaxis() -> SetTitleSize(0.05);
  h_mZ_DA[category] -> GetYaxis() -> SetTitleSize(0.071);
  h_mZ_DA[category] -> GetYaxis() -> SetTitleOffset(1.22);  
  h_mZ_DA[category] -> GetYaxis() -> SetTitle(axisTitle);
  
  h_mZ_DA[category] -> SetLineWidth(1);
  h_mZ_DA[category] -> SetLineColor(kBlack);
  h_mZ_DA[category] -> SetMarkerColor(kBlack);
  h_mZ_DA[category] -> SetMarkerStyle(20);
  h_mZ_DA[category] -> SetMarkerSize(0.7);
  gPad->Update();
  
  
  // fit the plots
  TF1* f_mee_voigtian;
  TF1* f_mee_crystalBall;
  FitZPeak(doVoigtianFit,doCrystalBallFit,&f_mee_voigtian,&f_mee_crystalBall,"f_mee",h_mZ_MC[category],EBEE);
  double peakV,sigmaOL,peakCB,sigmaCB;
  TArrow* line;
  TArrow* line2;
  GetOLSigma(peakV,sigmaOL,f_mee_voigtian,0.001,&line,&line2);
  peakCB = 90.1876+f_mee_crystalBall -> GetParameter(3);
  sigmaCB = f_mee_crystalBall -> GetParameter(4);
  
  
  char stringa[80];
  sprintf(stringa,"peak_{V} = %2.3f",peakV);
  TLatex* latex1 = new TLatex(0.20,0.85,stringa);
  latex1->SetNDC();
  latex1->SetTextSize(0.057);
  latex1->SetTextFont(42);
  latex1->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
  latex1->SetLineWidth(2);
  
  sprintf(stringa,"#sigma_{OL}/peak = %4.2f%s",sigmaOL/peakV*100.,"\%");
  TLatex* latex2 = new TLatex(0.20,0.77,stringa);
  latex2->SetNDC();
  latex2->SetTextSize(0.057);
  latex2->SetTextFont(42);
  latex2->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
  latex2->SetLineWidth(2);
  
  sprintf(stringa,"peak_{CB} = %2.3f",peakCB);
  TLatex* latex3 = new TLatex(0.20,0.67,stringa);
  latex3->SetNDC();
  latex3->SetTextSize(0.057);
  latex3->SetTextFont(42);
  latex3->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
  latex3->SetLineWidth(2);
  
  sprintf(stringa,"#sigma_{CB}/peak = %4.2f%s",sigmaCB/peakCB*100.,"\%");
  TLatex* latex4 = new TLatex(0.20,0.59,stringa);
  latex4->SetNDC();
  latex4->SetTextSize(0.057);
  latex4->SetTextFont(42);
  latex4->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
  latex4->SetLineWidth(2);
  
  
  
  TF1* f_c_mee_voigtian;
  TF1* f_c_mee_crystalBall;
  FitZPeak(doVoigtianFit,doCrystalBallFit,&f_c_mee_voigtian,&f_c_mee_crystalBall,"f_c_mee",h_mZ_DA[category],EBEE);
  
  double peakV_c,sigmaOL_c,peakCB_c,sigmaCB_c;
  TArrow* line_c;
  TArrow* line2_c;
  GetOLSigma(peakV_c,sigmaOL_c,f_c_mee_voigtian,0.001,&line_c,&line2_c);
  peakCB_c = 90.1876+f_c_mee_crystalBall -> GetParameter(3);
  sigmaCB_c = f_c_mee_crystalBall -> GetParameter(4);
  
  char stringa_c[80];
  sprintf(stringa_c,"peak_{V} = %2.3f",peakV_c);
  TLatex* latex1_c = new TLatex(0.62,0.85,stringa_c);
  latex1_c->SetNDC();
  latex1_c->SetTextSize(0.057);
  latex1_c->SetTextFont(42);
  latex1_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
  latex1_c->SetLineWidth(2);
  
  sprintf(stringa_c,"#sigma_{OL}/peak = %4.2f%s",sigmaOL_c/peakV_c*100.,"\%");
  TLatex* latex2_c = new TLatex(0.62,0.77,stringa_c);
  latex2_c->SetNDC();
  latex2_c->SetTextSize(0.057);
  latex2_c->SetTextFont(42);
  latex2_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
  latex2_c->SetLineWidth(2);
  
  sprintf(stringa_c,"peak_{CB} = %2.3f",peakCB_c);
  TLatex* latex3_c = new TLatex(0.62,0.67,stringa_c);
  latex3_c->SetNDC();
  latex3_c->SetTextSize(0.057);
  latex3_c->SetTextFont(42);
  latex3_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
  latex3_c->SetLineWidth(2);
  
  sprintf(stringa_c,"#sigma_{CB}/peak = %4.2f%s",sigmaCB_c/peakCB_c*100.,"\%");
  TLatex* latex4_c = new TLatex(0.62,0.59,stringa_c);
  latex4_c->SetNDC();
  latex4_c->SetTextSize(0.057);
  latex4_c->SetTextFont(42);
  latex4_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
  latex4_c->SetLineWidth(2);
  
  
  // print the plots
  p1 -> cd();
  
  float maximum = h_mZ_MC[category] -> GetMaximum();
  if( h_mZ_DA[category]->GetMaximum() > maximum) maximum = h_mZ_DA[category]->GetMaximum();
  
  h_mZ_MC[category] -> SetMinimum(0.);
  h_mZ_MC[category] -> SetMaximum(1.05*maximum);
  h_mZ_MC[category] -> Draw("hist");
  h_mZ_MC[category] -> Draw("P,same");
  h_mZ_DA[category] -> Draw("P,sames");
  
  f_mee_crystalBall -> Draw("same");
  f_mee_voigtian    -> Draw("same");
  f_c_mee_crystalBall -> Draw("same");
  f_c_mee_voigtian    -> Draw("same");
  
  line -> Draw("same");
  line_c -> Draw("same");
  
  line2 -> Draw("same");
  line2_c -> Draw("same");
  
  latex1   -> Draw("same");
  latex2   -> Draw("same");
  latex3   -> Draw("same");
  latex4   -> Draw("same");
  latex1_c -> Draw("same");
  latex2_c -> Draw("same");
  latex3_c -> Draw("same");
  latex4_c -> Draw("same");
  
  
  
  p2 -> cd();
  
  TArrow* line3;
  TH1F* ratio_MC = ratioHisto(h_mZ_MC[category],h_mZ_MC[category],peakV,&line3);
  TArrow* line4;
  TH1F* ratio_DA = ratioHisto(h_mZ_DA[category],h_mZ_MC[category],peakV_c,&line4);
  
  ratio_MC -> SetFillStyle(0);
  ratio_MC -> Draw();
  ratio_DA -> Draw("same");
  
  line3 -> Draw("same");
  line4 -> Draw("same");
  
  gPad -> Update();
  
  
  
  c -> Print((outFileName+"."+extension).c_str(),extension.c_str());
  delete c;
  
  outFile -> cd();
  
  h_mZ_MC[category] -> Write();
  h_mZ_DA[category] -> Write();
  
}



TH1F* ratioHisto(TH1F* h_num, TH1F* h_den, const double& xMax, TArrow** line)
{
  TH1F* h_ratio = (TH1F*)( h_num->Clone() );
  
  for(int bin = 1; bin <= h_num->GetNbinsX(); ++bin)
  {
    double int_num = h_num -> Integral(1,bin);
    double int_den = h_den -> Integral(1,bin);
    h_ratio -> SetBinContent(bin,int_num/int_den);
    h_ratio -> SetBinError(bin,0.);
  }
  
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
  h_ratio -> GetYaxis() -> SetTitle("cml ratio");
  
  if( (xMax != -999.) && (line != NULL) )
  {
    (*line) = new TArrow(xMax,yMin,xMax,yMax);
    (*line) -> SetLineWidth(2);
    (*line) -> SetLineColor(h_ratio->GetLineColor());
  }
  
  return h_ratio;
}
