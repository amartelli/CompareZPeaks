#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "geometryUtils.h"
#include "ZPeakFitUtils.h"
#include "PUReweighting.h"
#include "GetCategories.h"
#include "GetScaleCorrection.h"
#include "GetExtraSmearing.h"
#include "ParseIJazZFile.h"
#include "ScaleEstimators.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
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
#include "TGraphErrors.h"
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
               const bool& doVoigtianFit, const bool& doCrystalBallFit, const bool& doEffectiveSigma,
               const std::string& EBEE, const std::string& outFileName);

TH1F* ratioHisto(TH1F* h_num, TH1F* h_den, const double& xMax = -999., TArrow** line = NULL);



//-----------------
// global variables

double meeMin = 75.;
double meeMax = 105.;

std::map<std::string,TH1F*> h_mZ_MC;
std::map<std::string,TH1F*> h_mZ_DA;

double meeFineMin = 60.;
double meeFineMax = 120.;
double precision = 0.001;

std::map<std::string,TH1F*> h_mZFine_MC;
std::map<std::string,TH1F*> h_mZFine_DA;

std::map<std::string,TH1F*> h_mZRes_MC;

TFile* outFile;
std::string extension = "pdf";






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
  
  std::string dataLabel  = gConfigParser -> readStringOption("Options::dataLabel");
  
  bool useShervinNtuple = gConfigParser -> readBoolOption("Options::useShervinNtuple");
  
  int maxEntries    = gConfigParser -> readIntOption("Options::maxEntries");
  std::string MCGen = gConfigParser -> readStringOption("Options::MCGen");
  bool runDepFlag   = gConfigParser -> readBoolOption("Options::runDepFlag");
  int runMin        = gConfigParser -> readIntOption("Options::runMin");
  int runMax        = gConfigParser -> readIntOption("Options::runMax");
  
  bool diagonalCatOnly = gConfigParser -> readBoolOption("Options::diagonalCatOnly");
  
  std::string eleIDSelection = gConfigParser -> readStringOption("Options::eleIDSelection");
  int eleIDBit = 1;
  if( eleIDSelection == "loose"  ) eleIDBit = 2;
  if( eleIDSelection == "medium" ) eleIDBit = 6;
  if( eleIDSelection == "tight"  ) eleIDBit = 14;
  if( eleIDSelection == "WP90PU" ) eleIDBit = 16;
  if( eleIDSelection == "WP80PU" ) eleIDBit = 48;
    
  bool applyPUWeight        = gConfigParser -> readBoolOption("Options::applyPUWeight");
  bool applyEnergyScaleCorr = gConfigParser -> readBoolOption("Options::applyEnergyScaleCorr");
  bool applyEnergySmearing  = gConfigParser -> readBoolOption("Options::applyEnergySmearing");
  
  std::string enCorrType          = gConfigParser -> readStringOption("Options::enCorrType");
  std::string energyScaleCorrType = gConfigParser -> readStringOption("Options::energyScaleCorrType");
  std::string energySmearingType  = gConfigParser -> readStringOption("Options::energySmearingType");
  
  std::string runRangeFile        = gConfigParser -> readStringOption("Options::runRangeFile");
  std::string ShervinScaleFile    = gConfigParser -> readStringOption("Options::ShervinScaleFile");
  std::string ShervinSmearingFile = gConfigParser -> readStringOption("Options::ShervinSmearingFile");
  std::string IJazZGlobalFolder   = gConfigParser -> readStringOption("Options::IJazZGlobalFolder");
  std::string IJazZRunDepFolder   = gConfigParser -> readStringOption("Options::IJazZRunDepFolder");
  
  bool doVoigtianFit    = gConfigParser -> readBoolOption("Options::doVoigtianFit");
  bool doCrystalBallFit = gConfigParser -> readBoolOption("Options::doCrystalBallFit");
  bool doEffectiveSigma = gConfigParser -> readBoolOption("Options::doEffectiveSigma");
  
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
  std::cout << ">>> Get trees" << std::endl;  
  std::cout << std::endl;
  
  std::string treeName;
  if( !useShervinNtuple ) treeName = "simpleNtupleEoverP/SimpleNtupleEoverP";
  else                    treeName = "selected";
  
  TChain* ntu_MC = new TChain(treeName.c_str());
  FillChain(ntu_MC,inputFilesMC);
  std::cout << ">>>   MC: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  
  TChain* ntu_DA = new TChain(treeName.c_str());
  FillChain(ntu_DA,inputFilesDA);
  std::cout << ">>> DATA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;
  
  if( ntu_MC->GetEntries() == 0 || ntu_DA->GetEntries() == 0 )
  {
    std::cout << ">>> compareZPeaks::Error: at least one file is empty" << std::endl; 
    return -1;
  }
  
  
  
  //------------------------
  // Define branch addresses
  std::cout << std::endl;
  std::cout << ">>> Define branch addresses" << std::endl;
  
  float scEta[2];
  float scE[2];
  float scERaw[2];
  float scEReg[2];
  float scETrue[2];
  float R9[2];
  int eleID[2];
  
  int runId;
  int isZ;
  bool HLTfire;
  int timeStampHigh;
  float mZ, mZTrue;
  int isEB1,isEB2;
  float scEta1,scEta2;
  float scERaw1,scERaw2;
  float scE1,scE2;
  float scEReg1,scEReg2;
  float scETrue1,scETrue2;
  float E3x31,E3x32;
  float R91,R92;
  int seedIeta1,seedIeta2;
  int seedIphi1,seedIphi2;
  int seedIx1,seedIx2;
  int seedIy1,seedIy2;
  int seedIz1,seedIz2;
  int eleID1, eleID2;
  int nPU;
  
  if( !useShervinNtuple )
  {
    HLTfire = true;
    
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
    ntu_MC -> SetBranchStatus("PUit_TrueNumInteractions",1); ntu_MC -> SetBranchAddress("PUit_TrueNumInteractions",&nPU);  
  }
  else
  {
    isZ = 1;
    
    ntu_MC -> SetBranchStatus("*",0);
    ntu_MC -> SetBranchStatus("HLTfire",       1);   ntu_MC -> SetBranchAddress("HLTfire",&HLTfire);
    ntu_MC -> SetBranchStatus("runNumber",     1);   ntu_MC -> SetBranchAddress("runNumber",&runId);
    ntu_MC -> SetBranchStatus("nPU",           1);   ntu_MC -> SetBranchAddress("nPU",&nPU);
    ntu_MC -> SetBranchStatus("R9Ele",         1);   ntu_MC -> SetBranchAddress("R9Ele",R9);
    ntu_MC -> SetBranchStatus("etaSCEle",      1);   ntu_MC -> SetBranchAddress("etaSCEle",scEta);
    ntu_MC -> SetBranchStatus("rawEnergySCEle",1);   ntu_MC -> SetBranchAddress("rawEnergySCEle",scERaw);
    ntu_MC -> SetBranchStatus("energyMCEle",   1);   ntu_MC -> SetBranchAddress("energyMCEle",scETrue);
    ntu_MC -> SetBranchStatus("energySCEle",   1);   ntu_MC -> SetBranchAddress("energySCEle",scE);
    if( enCorrType == "stdSC" )
    {
      ntu_MC -> SetBranchStatus("energySCEle",1);   ntu_MC -> SetBranchAddress("energySCEle",scEReg);
    }
    if( enCorrType == "eleTunedRegV3" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorr_ele",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorr_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV3" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorr_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorr_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV4" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV4_ele",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV4_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV4" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV4_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV4_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV5" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV5_ele",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV5_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV5" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV5_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV5_pho",scEReg);
    }
    ntu_MC -> SetBranchStatus("invMass_SC",1);   ntu_MC -> SetBranchAddress("invMass_SC",&mZ);
    ntu_MC -> SetBranchStatus("eleID",     1);   ntu_MC -> SetBranchAddress("eleID",eleID);
    
    ntu_DA -> SetBranchStatus("*",0);
    ntu_DA -> SetBranchStatus("HLTfire",       1);   ntu_DA -> SetBranchAddress("HLTfire",&HLTfire);
    ntu_DA -> SetBranchStatus("runNumber",     1);   ntu_DA -> SetBranchAddress("runNumber",&runId);
    ntu_DA -> SetBranchStatus("R9Ele",         1);   ntu_DA -> SetBranchAddress("R9Ele",R9);
    ntu_DA -> SetBranchStatus("etaSCEle",      1);   ntu_DA -> SetBranchAddress("etaSCEle",scEta);
    ntu_DA -> SetBranchStatus("rawEnergySCEle",1);   ntu_DA -> SetBranchAddress("rawEnergySCEle",scERaw);
    ntu_DA -> SetBranchStatus("energySCEle",   1);   ntu_DA -> SetBranchAddress("energySCEle",scE);
    if( enCorrType == "stdSC" )
    {
      ntu_DA -> SetBranchStatus("energySCEle",1);   ntu_DA -> SetBranchAddress("energySCEle",scEReg);
    }
    if( enCorrType == "eleTunedRegV3" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorr_ele",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorr_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV3" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorr_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorr_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV4" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV4_ele",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV4_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV4" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV4_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV4_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV5" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV5_ele",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV5_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV5" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV5_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV5_pho",scEReg);
    }
    ntu_DA -> SetBranchStatus("invMass_SC",1);   ntu_DA -> SetBranchAddress("invMass_SC", &mZ);
    ntu_DA -> SetBranchStatus("eleID",     1);   ntu_DA -> SetBranchAddress("eleID",eleID);
  }
  
  
  
  //------------------
  // Define categories
  std::cout << std::endl;
  std::cout << ">>> Define categories" << std::endl;
    
  std::vector<std::string> categories;
  std::vector<std::string> categoryLabels;
  
  categories.push_back("EB-EB");
  categoryLabels.push_back("EB-EB");
  categories.push_back("EB-EB_hR9");
  categoryLabels.push_back("EB-EB   R9>0.94");
  categories.push_back("EB-EB_lR9");
  categoryLabels.push_back("EB-EB   R9<0.94");
  categories.push_back("EB-EB_eta>1");
  categoryLabels.push_back("EB-EB   |#eta|>1");
  categories.push_back("EB-EB_eta<1");
  categoryLabels.push_back("EB-EB   |#eta|<1");
  
  categories.push_back("EE-EE");
  categoryLabels.push_back("EE-EE");
  categories.push_back("EE-EE_hR9");
  categoryLabels.push_back("EE-EE   R9>0.94");
  categories.push_back("EE-EE_lR9");
  categoryLabels.push_back("EE-EE   R9<0.94");
  categories.push_back("EE-EE_eta>2");
  categoryLabels.push_back("EE-EE   |#eta|>2");
  categories.push_back("EE-EE_eta<2");
  categoryLabels.push_back("EE-EE   |#eta|<2");
  
  std::vector<float> etaBins;
  etaBins.push_back(0.0);
  etaBins.push_back(1.0);
  etaBins.push_back(1.5);
  etaBins.push_back(2.0);
  etaBins.push_back(2.5);
  int nEtaBins = int(etaBins.size()) - 1;
  
  std::vector<float> R9Bins;
  R9Bins.push_back(0.00);
  R9Bins.push_back(0.94);
  R9Bins.push_back(1.00);
  int nR9Bins = int(R9Bins.size()) - 1;
  
  std::vector<std::string> singleCats;
  std::vector<std::string> singleCatLabels;
  for(int etaBin = 0; etaBin < nEtaBins; ++etaBin)
    for(int R9Bin = 0; R9Bin < nR9Bins; ++R9Bin)
    {
      singleCats.push_back(Form("eta%1.1f-%1.1f_R9%1.2f-%1.2f",etaBins.at(etaBin),etaBins.at(etaBin+1),R9Bins.at(R9Bin),R9Bins.at(R9Bin+1)));
      singleCatLabels.push_back(Form("%1.1f<|#eta|<%1.1f %1.2f<R9<%1.2f",etaBins.at(etaBin),etaBins.at(etaBin+1),R9Bins.at(R9Bin),R9Bins.at(R9Bin+1)));
    }
  
  for(unsigned int i = 0; i < singleCats.size(); ++i)
    for(unsigned int j = i; j < singleCats.size(); ++j)
    {
      if( diagonalCatOnly && (j != i ) ) continue;
      
      categories.push_back(singleCats.at(i) + "__" + singleCats.at(j));
      categoryLabels.push_back(Form("%s %s",singleCatLabels.at(i).c_str(),singleCatLabels.at(j).c_str()));
      std::cout << ">>>>>> " << singleCats.at(i) + "__" + singleCats.at(j) << std::endl;
    }
  
  
  
  //------------------
  // Define histograms
  std::cout << std::endl;
  std::cout << ">>> Define histograms" << std::endl;
    
  for(unsigned int i = 0; i < categories.size(); ++i)
  {
    std::string category = categories.at(i);
    
    std::string histoName = "h_mZ_MC_"+category;
    h_mZ_MC[category] = new TH1F(histoName.c_str(),"",100,meeMin,meeMax);
    h_mZ_MC[category] -> Sumw2();
    
    histoName = "h_mZFine_MC_"+category;
    h_mZFine_MC[category] = new TH1F(histoName.c_str(),"",int((meeFineMax-meeFineMin)/precision),meeFineMin,meeFineMax);
    h_mZFine_MC[category] -> Sumw2();
    
    histoName = "h_mZRes_MC_"+category;
    h_mZRes_MC[category] = new TH1F(histoName.c_str(),"",10000,0.,2.);
    h_mZRes_MC[category] -> Sumw2();
    
    histoName = "h_mZ_DA_"+category;
    h_mZ_DA[category] = new TH1F(histoName.c_str(),"",100,meeMin,meeMax);
    h_mZ_DA[category] -> Sumw2();
    
    histoName = "h_mZFine_DA_"+category;
    h_mZFine_DA[category] = new TH1F(histoName.c_str(),"",int((meeFineMax-meeFineMin)/precision),meeFineMin,meeFineMax);
    h_mZFine_DA[category] -> Sumw2();
  }
  
  
  
  //-----------------------------
  // Setup data scale corrections
  std::cout << std::endl;
  std::cout << ">>> Setup data scale corrections" << std::endl;
  
  ScaleCorrector* myScaleCorrector = new ScaleCorrector(runRangeFile);
  
  if( applyEnergyScaleCorr )
  {
    if( energyScaleCorrType == "shervin" ) myScaleCorrector -> SetShervinRunDepScaleMap(ShervinScaleFile);
    if( energyScaleCorrType == "fabrice" ) myScaleCorrector -> SetIJazZGlobalScaleHisto(IJazZGlobalFolder);
    if( energyScaleCorrType == "fabrice" ) myScaleCorrector -> SetIJazZRunDepScaleHistoMap(IJazZRunDepFolder);
  }
  
  
  
  //-----------------------
  // Setup MC extrasmearing
  std::cout << std::endl;
  std::cout << ">>> Setup MC extrasmearing" << std::endl;
  
  Smearer* mySmearer = new Smearer();
  
  if( applyEnergySmearing )
  {
    if( energyScaleCorrType == "shervin" ) mySmearer -> SetShervinExtraSmearingMap(ShervinSmearingFile);
    if( energyScaleCorrType == "fabrice" ) mySmearer -> SetIJazZExtraSmearingHisto(IJazZGlobalFolder);
  }
  
  
  
  //--------------------------
  // pileup reweighting for MC
  std::cout << std::endl;
  std::cout << ">>> Setup MC pileup reweighting" << std::endl;
  
  std::map<std::string, TH1F*>* PUWeights = ReadPUWeights(MCGen,runDepFlag,runMin,runMax);
  
  
  
  
  
  
  //------------------
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
    R91 = E3x31/scERaw1;
    R92 = E3x32/scERaw2;
    
    if( useShervinNtuple )
    {
      R91 = R9[0];
      R92 = R9[1];
      scEta1 = scEta[0];
      scEta2 = scEta[1];
      scEReg1 = scEReg[0];
      scEReg2 = scEReg[1];
      scE1 = scE[0];
      scE2 = scE[1];
      eleID1 = eleID[0];
      eleID2 = eleID[1];
      scETrue1 = scETrue[0];
      scETrue2 = scETrue[1];
    }
    else
    {
      eleID1 = 7;
      eleID2 = 7;
    }
    
    
    // selections
    if( isZ != 1 ) continue;
    if( !HLTfire ) continue;
    if( fabs(scEta1) >= 2.5000 || fabs(scEta2) >= 2.5000  ) continue;
    if( fabs(scEta1) >  1.4442 && fabs(scEta1) <  1.5660 ) continue;
    if( fabs(scEta2) >  1.4442 && fabs(scEta2) <  1.5660 ) continue;
    if( R91 < 0.0 || R91 >= 1.0 ) continue;
    if( R92 < 0.0 || R92 >= 1.0 ) continue;
    if( ((eleID1 & eleIDBit) != eleIDBit) || ((eleID2 & eleIDBit) != eleIDBit) ) continue;
    
    
    if( applyEnergySmearing )
    {
      float energySmearing1 = gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta1,R91,dataLabel,energySmearingType));
      float energySmearing2 = gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta2,R92,dataLabel,energySmearingType));
      scEReg1 *= energySmearing1;
      scEReg2 *= energySmearing2;
    }
    
    float ww = 1.;
    if( applyPUWeight )
    {
      std::string periodLabel = getPeriodLabel(runId,runDepFlag,runMin,runMax);
      
      int ibin = (*PUWeights)[periodLabel] -> FindBin( nPU );
      if( ibin <= 1 ) ibin = 1;
      if( ibin >= (*PUWeights)[periodLabel]->GetNbinsX() ) ibin = (*PUWeights)[periodLabel]->GetNbinsX();
      ww *= (*PUWeights)[periodLabel]->GetBinContent(ibin);
    }
    
    
    // use regression energy
    mZTrue = mZ * sqrt( scETrue1/scE1 * scETrue2/scE2 ); 
    mZ *= sqrt( scEReg1/scE1 * scEReg2/scE2 );
            
    // fill EB histograms
    if( (fabs(scEta1) < 1.5) && (fabs(scEta2) < 1.5) )
    {
      h_mZ_MC["EB-EB"] -> Fill( mZ,ww );
      h_mZFine_MC["EB-EB"] -> Fill( mZ,ww );
      h_mZRes_MC["EB-EB"] -> Fill( mZ/mZTrue,ww );
      
      if( (R91 > 0.94) && (R92 > 0.94) )
      {
        h_mZ_MC["EB-EB_hR9"] -> Fill( mZ,ww );
        h_mZFine_MC["EB-EB_hR9"] -> Fill( mZ,ww );
        h_mZRes_MC["EB-EB_hR9"] -> Fill( mZ/mZTrue,ww );
      }
      
      if( (R91 < 0.94) && (R92 < 0.94) )
      {
        h_mZ_MC["EB-EB_lR9"] -> Fill( mZ,ww );
        h_mZFine_MC["EB-EB_lR9"] -> Fill( mZ,ww );
        h_mZRes_MC["EB-EB_lR9"] -> Fill( mZ/mZTrue,ww );
      }
      
      if( (fabs(scEta1) > 1.) && (fabs(scEta2) > 1.) )
      {
        h_mZ_MC["EB-EB_eta>1"] -> Fill( mZ,ww );
        h_mZFine_MC["EB-EB_eta>1"] -> Fill( mZ,ww );
        h_mZRes_MC["EB-EB_eta>1"] -> Fill( mZ/mZTrue,ww );
      }
      
      if( (fabs(scEta1) < 1.) && (fabs(scEta2) < 1.) )
      {
        h_mZ_MC["EB-EB_eta<1"] -> Fill( mZ,ww );
        h_mZFine_MC["EB-EB_eta<1"] -> Fill( mZ,ww );
        h_mZRes_MC["EB-EB_eta<1"] -> Fill( mZ/mZTrue,ww );
      }
    }
    
    
    // fill EE histograms
    if( (fabs(scEta1) > 1.5) && (fabs(scEta2) > 1.5) )
    {
      h_mZ_MC["EE-EE"] -> Fill( mZ,ww );
      h_mZFine_MC["EE-EE"] -> Fill( mZ,ww );
      h_mZRes_MC["EE-EE"] -> Fill( mZ/mZTrue,ww );
      
      if( (R91 > 0.94) && (R92 > 0.94) )
      {
        h_mZ_MC["EE-EE_hR9"] -> Fill( mZ,ww );
        h_mZFine_MC["EE-EE_hR9"] -> Fill( mZ,ww );
        h_mZRes_MC["EE-EE_hR9"] -> Fill( mZ/mZTrue,ww );
      }
      
      if( (R91 < 0.94) && (R92 < 0.94) )
      {
        h_mZ_MC["EE-EE_lR9"] -> Fill( mZ,ww );
        h_mZFine_MC["EE-EE_lR9"] -> Fill( mZ,ww );
        h_mZRes_MC["EE-EE_lR9"] -> Fill( mZ/mZTrue,ww );
      }
      
      if( (fabs(scEta1) > 2.) && (fabs(scEta2) > 2.) )
      {
        h_mZ_MC["EE-EE_eta>2"] -> Fill( mZ,ww );
        h_mZFine_MC["EE-EE_eta>2"] -> Fill( mZ,ww );
        h_mZRes_MC["EE-EE_eta>2"] -> Fill( mZ/mZTrue,ww );
      }
      
      if( (fabs(scEta1) < 2.) && (fabs(scEta2) < 2.) )
      {
        h_mZ_MC["EE-EE_eta<2"] -> Fill( mZ,ww );
        h_mZFine_MC["EE-EE_eta<2"] -> Fill( mZ,ww );
        h_mZRes_MC["EE-EE_eta<2"] -> Fill( mZ/mZTrue,ww );
      }
    }
    
    
    // fill all independent categories
    std::string catLabel = GetEtaR9CatLabel(fabs(scEta1),R91,fabs(scEta2),R92,etaBins,R9Bins,categories);
    if( catLabel != "undefined__undefined" )
    {
      h_mZ_MC[catLabel] -> Fill( mZ,ww );
      h_mZFine_MC[catLabel] -> Fill( mZ,ww );
      h_mZRes_MC[catLabel] -> Fill( mZ/mZTrue,ww );
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
    
    if( useShervinNtuple )
    {
      R91 = R9[0];
      R92 = R9[1];
      scEta1 = scEta[0];
      scEta2 = scEta[1];
      scEReg1 = scEReg[0];
      scEReg2 = scEReg[1];
      scE1 = scE[0];
      scE2 = scE[1];
      eleID1 = eleID[0];
      eleID2 = eleID[1];
    }
    else
    {
      eleID1 = 7;
      eleID2 = 7;
    }
    
    
    // selections
    if( isZ != 1 ) continue;
    if( !HLTfire ) continue;
    if( fabs(scEta1) >= 2.5000 || fabs(scEta2) >= 2.5000  ) continue;
    if( fabs(scEta1) >  1.4442 && fabs(scEta1) <  1.5660 ) continue;
    if( fabs(scEta2) >  1.4442 && fabs(scEta2) <  1.5660 ) continue;
    if( R91 < 0.0 || R91 >= 1.0 ) continue;
    if( R92 < 0.0 || R92 >= 1.0 ) continue;
    if( ((eleID1 & eleIDBit) != eleIDBit) || ((eleID2 & eleIDBit) != eleIDBit) ) continue;    
    
    
    if( applyEnergyScaleCorr )
    {
      scEReg1 *= myScaleCorrector->GetScaleCorrection(scEta1,R91,runId,dataLabel,energyScaleCorrType);
      scEReg2 *= myScaleCorrector->GetScaleCorrection(scEta2,R92,runId,dataLabel,energyScaleCorrType);
    }
    
    
    // use regression energy
    mZ *= sqrt( scEReg1/scE1 * scEReg2/scE2 );
    
    
    // fill EB histograms
    if( (fabs(scEta1) < 1.5) && (fabs(scEta2) < 1.5) )
    {
      h_mZ_DA["EB-EB"]  ->  Fill( mZ );
      h_mZFine_DA["EB-EB"]  ->  Fill( mZ );
      
      if( (R91 > 0.94) && (R92 > 0.94) )
      {
        h_mZ_DA["EB-EB_hR9"]  ->  Fill( mZ );
        h_mZFine_DA["EB-EB_hR9"]  ->  Fill( mZ );
      }
      
      if( (R91 < 0.94) && (R92 < 0.94) )
      {
        h_mZ_DA["EB-EB_lR9"]  ->  Fill( mZ );
        h_mZFine_DA["EB-EB_lR9"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 1.) && (fabs(scEta2) > 1.) )
      {
        h_mZ_DA["EB-EB_eta>1"]  ->  Fill( mZ );
        h_mZFine_DA["EB-EB_eta>1"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 1.) && (fabs(scEta2) < 1.) )
      {
        h_mZ_DA["EB-EB_eta<1"]  ->  Fill( mZ );
        h_mZFine_DA["EB-EB_eta<1"]  ->  Fill( mZ );
      }
    }
    
    
    // fill EE histograms
    if( (fabs(scEta1) > 1.5) && (fabs(scEta2) > 1.5) )
    {
      h_mZ_DA["EE-EE"]  ->  Fill( mZ );
      h_mZFine_DA["EE-EE"]  ->  Fill( mZ );
      
      if( (R91 > 0.94) && (R92 > 0.94) )
      {
        h_mZ_DA["EE-EE_hR9"]  ->  Fill( mZ );
        h_mZFine_DA["EE-EE_hR9"]  ->  Fill( mZ );
      }
      
      if( (R91 < 0.94) && (R92 < 0.94) )
      {
        h_mZ_DA["EE-EE_lR9"]  ->  Fill( mZ );
        h_mZFine_DA["EE-EE_lR9"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) > 2.) && (fabs(scEta2) > 2.) )
      {
        h_mZ_DA["EE-EE_eta>2"]  ->  Fill( mZ );
        h_mZFine_DA["EE-EE_eta>2"]  ->  Fill( mZ );
      }
      
      if( (fabs(scEta1) < 2.) && (fabs(scEta2) < 2.) )
      {
        h_mZ_DA["EE-EE_eta<2"]  ->  Fill( mZ );
        h_mZFine_DA["EE-EE_eta<2"]  ->  Fill( mZ );
      }
    }
    
    
    // fill all independent categories
    std::string catLabel = GetEtaR9CatLabel(fabs(scEta1),R91,fabs(scEta2),R92,etaBins,R9Bins,categories);
    if( catLabel != "undefined__undefined" )
    {
      h_mZ_DA[catLabel] -> Fill( mZ );
      h_mZFine_DA[catLabel] -> Fill( mZ );
    }
  }
  std::cout << std::endl;
  
  
  
  //---------
  // Drawings
  std::cout << ">>> Drawings" << std::endl;
  
  std::string folderName = outFilePath + "/" + dataLabel + "/";
  gSystem -> mkdir(folderName.c_str());
  if( energySmearingType == "fabrice" ) folderName += IJazZGlobalFolder + "/";
  gSystem -> mkdir(folderName.c_str());
  
  std::string outFileName = folderName + "/compareZPeaks";
  
  outFileName += "__" + enCorrType;
  outFileName += "__MC_" + MCGen;
  if( applyPUWeight )        outFileName += "__MC_PUweight";
  if( applyEnergySmearing )  outFileName += "__MC_energySmearing_" + energySmearingType;
  if( applyEnergyScaleCorr ) outFileName += "__DA_energyScaleCorr_" + energyScaleCorrType;
  
  outFile = new TFile((outFileName+".root").c_str(),"RECREATE");
  
  
  TCanvas* dummy = new TCanvas("dummy","",0,0,800,600);
  dummy -> Print((outFileName+"."+extension+"[").c_str(),extension.c_str());
  
  for(unsigned int cat = 0; cat < categories.size(); ++cat)
  {
    std::string category = categories.at(cat);
    DrawZPeak(categories.at(cat),categoryLabels.at(cat),1,doVoigtianFit,doCrystalBallFit,doEffectiveSigma,"EBEE",outFileName);
  }
  
  outFile -> Close();
  dummy -> Print((outFileName+"."+extension+"]").c_str(),extension.c_str());}







void DrawZPeak(const std::string& category, const std::string& label, const int& rebin,
               const bool& doVoigtianFit, const bool& doCrystalBallFit, const bool& doEffectiveSigma,
               const std::string& EBEE, const std::string& outFileName)
{
  TCanvas* c = new TCanvas(("c_mZ_"+category).c_str(),("mZ - "+category).c_str(),0,0,800,600);
  c -> cd();
  
  TPad* p1 = new TPad("p1","p1",0.00, 0.30, 0.75, 1.00);
  TPad* p2 = new TPad("p2","p2",0.00, 0.00, 0.75, 0.30);
  TPad* p3 = new TPad("p3","p3",0.75, 0.30, 1.00, 1.00);
  TPad* p4 = new TPad("p4","p4",0.75, 0.00, 1.00, 0.30);
  p1 -> SetLeftMargin(0.16);
  p1 -> SetRightMargin(0.01);
  p2 -> SetLeftMargin(0.16);
  p2 -> SetRightMargin(0.01);
  p4 -> SetLeftMargin(0.06);
  p4 -> SetRightMargin(0.23);
  p1 -> SetTopMargin(0.08);
  p1 -> SetBottomMargin(0.02);
  p2 -> SetTopMargin(0.04);
  p2 -> SetBottomMargin(0.4);
  p4 -> SetTopMargin(0.04);
  p4 -> SetBottomMargin(0.4);
  
  p1 -> Draw();
  p2 -> Draw();
  p3 -> Draw();
  p4 -> Draw();
  
  p1 -> cd();
  p1 -> SetGridx();
  
  
  //-----------
  // draw peaks
  std::cout << ">>> compareZPeaks::DrawZPeak::draw peaks for category " << category << std::endl;
  
  if( h_mZ_MC[category]->Integral() > 0 )
    h_mZ_MC[category] -> Scale( 1. * h_mZ_DA[category]->Integral() / h_mZ_MC[category]->Integral() );
  
  char axisTitle[50];
  h_mZ_MC[category] -> Rebin(rebin);
  sprintf(axisTitle,"events / %.2e GeV",h_mZ_MC[category]->GetBinWidth(1));
  //h_mZ_MC[category] -> GetXaxis() -> SetLabelSize(0.04);
  h_mZ_MC[category] -> GetXaxis() -> SetLabelSize(0.);
  h_mZ_MC[category] -> GetXaxis() -> SetLabelFont(42);
  //h_mZ_MC[category] -> GetXaxis() -> SetTitleSize(0.05);
  h_mZ_MC[category] -> GetXaxis() -> SetTitleSize(0.);
  h_mZ_MC[category] -> GetXaxis() -> SetTitleOffset(1.20);
  h_mZ_MC[category] -> GetXaxis() -> SetTitle(("m_{ee} (GeV)   "+label).c_str());
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
  sprintf(axisTitle,"events / %.2e GeV",h_mZ_DA[category]->GetBinWidth(1));
  //h_mZ_DA[category] -> GetXaxis() -> SetLabelSize(0.04);
  h_mZ_DA[category] -> GetXaxis() -> SetLabelSize(0.);
  h_mZ_DA[category] -> GetXaxis() -> SetLabelFont(42);
  //h_mZ_DA[category] -> GetXaxis() -> SetTitleSize(0.05);
  h_mZ_DA[category] -> GetXaxis() -> SetTitleSize(0.);
  h_mZ_DA[category] -> GetXaxis() -> SetTitleOffset(1.20);  
  h_mZ_DA[category] -> GetXaxis() -> SetTitle(("m_{ee} (GeV)   "+label).c_str());
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
  
  TLatex* latex1;
  TLatex* latex2;
  TLatex* latex3;
  TLatex* latex4;
  char stringa[80];
  if( doVoigtianFit )
  {
    sprintf(stringa,"peak_{V} = %2.3f",peakV);
    latex1 = new TLatex(0.20,0.85,stringa);
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
    latex1->SetLineWidth(2);
    
    sprintf(stringa,"#sigma_{OL}/p. = %4.2f%s",sigmaOL/peakV*100.,"\%");
    latex2 = new TLatex(0.20,0.77,stringa);
    latex2->SetNDC();
    latex2->SetTextSize(0.05);
    latex2->SetTextFont(42);
    latex2->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
    latex2->SetLineWidth(2);
  }
  if( doCrystalBallFit )
  {
    sprintf(stringa,"peak_{CB} = %2.3f",peakCB);
    latex3 = new TLatex(0.20,0.67,stringa);
    latex3->SetNDC();
    latex3->SetTextSize(0.05);
    latex3->SetTextFont(42);
    latex3->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
    latex3->SetLineWidth(2);
    
    sprintf(stringa,"#sigma_{CB}/p. = %4.2f%s",sigmaCB/peakCB*100.,"\%");
    latex4 = new TLatex(0.20,0.59,stringa);
    latex4->SetNDC();
    latex4->SetTextSize(0.05);
    latex4->SetTextFont(42);
    latex4->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
    latex4->SetLineWidth(2);
  }
  
  
  TF1* f_c_mee_voigtian;
  TF1* f_c_mee_crystalBall;
  FitZPeak(doVoigtianFit,doCrystalBallFit,&f_c_mee_voigtian,&f_c_mee_crystalBall,"f_c_mee",h_mZ_DA[category],EBEE);
  
  double peakV_c,sigmaOL_c,peakCB_c,sigmaCB_c;
  TArrow* line_c;
  TArrow* line2_c;
  GetOLSigma(peakV_c,sigmaOL_c,f_c_mee_voigtian,0.001,&line_c,&line2_c);
  peakCB_c = 90.1876+f_c_mee_crystalBall -> GetParameter(3);
  sigmaCB_c = f_c_mee_crystalBall -> GetParameter(4);

  TLatex* latex1_c;
  TLatex* latex2_c;
  TLatex* latex3_c;
  TLatex* latex4_c;  
  char stringa_c[80];
  if( doVoigtianFit )
  {
    sprintf(stringa_c,"peak_{V} = %2.3f",peakV_c);
    latex1_c = new TLatex(0.75,0.85,stringa_c);
    latex1_c->SetNDC();
    latex1_c->SetTextSize(0.05);
    latex1_c->SetTextFont(42);
    latex1_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
    latex1_c->SetLineWidth(2);
    
    sprintf(stringa_c,"#sigma_{OL}/p. = %4.2f%s",sigmaOL_c/peakV_c*100.,"\%");
    latex2_c = new TLatex(0.75,0.77,stringa_c);
    latex2_c->SetNDC();
    latex2_c->SetTextSize(0.05);
    latex2_c->SetTextFont(42);
    latex2_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
    latex2_c->SetLineWidth(2);
  }
  if( doCrystalBallFit )
  {
    sprintf(stringa_c,"peak_{CB} = %2.3f",peakCB_c);
    latex3_c = new TLatex(0.75,0.67,stringa_c);
    latex3_c->SetNDC();
    latex3_c->SetTextSize(0.05);
    latex3_c->SetTextFont(42);
    latex3_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
    latex3_c->SetLineWidth(2);
    
    sprintf(stringa_c,"#sigma_{CB}/p. = %4.2f%s",sigmaCB_c/peakCB_c*100.,"\%");
    latex4_c = new TLatex(0.75,0.59,stringa_c);
    latex4_c->SetNDC();
    latex4_c->SetTextSize(0.05);
    latex4_c->SetTextFont(42);
    latex4_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
    latex4_c->SetLineWidth(2);
  }
  
  
  // get the effective sigma
  TLatex* latex5;
  TLatex* latex5_c;
  TLatex* latex6;
  TLatex* latex6_c;
  
  if( doEffectiveSigma )
  {
    double mean,meanErr,min,max;
    double mean_c,meanErr_c,min_c,max_c;
    
    FindSmallestInterval(mean,  meanErr,  min,  max,  h_mZFine_MC[category],0.68);
    FindSmallestInterval(mean_c,meanErr_c,min_c,max_c,h_mZFine_DA[category],0.68);
    
    sprintf(stringa,"#sigma_{eff}^{68%s}/p. = %1.2f%s","\%",0.5*(max-min)/mean*100,"\%");
    latex5 = new TLatex(0.20,0.55,stringa);
    latex5->SetNDC();
    latex5->SetTextSize(0.05);
    latex5->SetTextFont(42);
    latex5->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
    latex5->SetLineWidth(2);
    
    sprintf(stringa_c,"#sigma_{eff}^{68%s}/p. = %1.2f%s","\%",0.5*(max_c-min_c)/mean_c*100,"\%");
    latex5_c = new TLatex(0.75,0.55,stringa_c);
    latex5_c->SetNDC();
    latex5_c->SetTextSize(0.05);
    latex5_c->SetTextFont(42);
    latex5_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
    latex5_c->SetLineWidth(2);
    
    FindSmallestInterval(mean,   meanErr,  min,  max,  h_mZFine_MC[category],0.38);
    FindSmallestInterval(mean_c,meanErr_c, min_c,max_c,h_mZFine_DA[category],0.38);
    
    sprintf(stringa,"#sigma_{eff}^{38%s}/p. = %1.2f%s","\%",0.5*(max-min)/mean*100,"\%");
    latex6 = new TLatex(0.20,0.47,stringa);
    latex6->SetNDC();
    latex6->SetTextSize(0.05);
    latex6->SetTextFont(42);
    latex6->SetTextColor(h_mZ_MC[category]->GetMarkerColor());
    latex6->SetLineWidth(2);
    
    sprintf(stringa_c,"#sigma_{eff}^{38%s}/p. = %1.2f%s","\%",0.5*(max_c-min_c)/mean_c*100,"\%");
    latex6_c = new TLatex(0.75,0.47,stringa_c);
    latex6_c->SetNDC();
    latex6_c->SetTextSize(0.05);
    latex6_c->SetTextFont(42);
    latex6_c->SetTextColor(h_mZ_DA[category]->GetMarkerColor());
    latex6_c->SetLineWidth(2);
  }
  
  
  // print the plots
  p1 -> cd();
  
  float maximum = h_mZ_MC[category] -> GetMaximum();
  if( h_mZ_DA[category]->GetMaximum() > maximum) maximum = h_mZ_DA[category]->GetMaximum();
  
  h_mZ_MC[category] -> SetMinimum(0.);
  h_mZ_MC[category] -> SetMaximum(1.05*maximum);
  h_mZ_MC[category] -> Draw("hist");
  h_mZ_MC[category] -> Draw("P,same");
  h_mZ_DA[category] -> Draw("P,sames");
  
  if(doCrystalBallFit) f_mee_crystalBall   -> Draw("same");
  if(doCrystalBallFit) f_c_mee_crystalBall -> Draw("same");
  if(doVoigtianFit)    f_mee_voigtian   -> Draw("same");
  if(doVoigtianFit)    f_c_mee_voigtian -> Draw("same");
  
  line -> Draw("same");
  line_c -> Draw("same");
  
  line2 -> Draw("same");
  line2_c -> Draw("same");
  
  if( doVoigtianFit )
  {
    latex1   -> Draw("same");
    latex2   -> Draw("same");
    latex1_c -> Draw("same");
    latex2_c -> Draw("same");
  }
  if( doCrystalBallFit )
  {
    latex3   -> Draw("same");
    latex4   -> Draw("same");
    latex3_c -> Draw("same");
    latex4_c -> Draw("same");
  }
  if( doEffectiveSigma )
  {
    latex5   -> Draw("same");
    latex5_c -> Draw("same");
    latex6   -> Draw("same");
    latex6_c -> Draw("same");
  }
  
  
  p2 -> cd();
  p2 -> SetGridx();
  
  double* residuals = new double[h_mZ_DA[category]->GetNbinsX()];
  double redChi2 = h_mZ_DA[category] -> Chi2Test(h_mZ_MC[category],"UW,CHI2/NDF",residuals);
  
  TH1F* h_residuals = (TH1F*)( h_mZ_DA[category]->Clone("h_residuals") );
  h_residuals -> Reset();
  
  TH1F* h_residualDistr = new TH1F("h_residualDistr","",20,-5,5);
  
  int ndf = 0;
  for(int bin = 1; bin <= h_mZ_DA[category]->GetNbinsX(); ++bin)
  {
    h_residuals -> SetBinContent(bin,residuals[bin-1]);
    
    if( !(h_mZ_DA[category]->GetBinContent(bin) == 0 && h_mZ_MC[category]->GetBinContent(bin) == 0) )
    {
      h_residualDistr -> Fill(residuals[bin-1]);
      ++ndf;
    }
  }
  
  h_residuals -> GetXaxis() -> SetLabelSize(0.13);
  h_residuals -> GetXaxis() -> SetLabelFont(42);
  h_residuals -> GetXaxis() -> SetTitleSize(0.17);
  h_residuals -> GetXaxis() -> SetTitleOffset(1.00);
  h_residuals -> GetXaxis() -> SetTitle("m_{ee}");
  h_residuals -> GetYaxis() -> SetLabelSize(0.13);
  h_residuals -> GetYaxis() -> SetLabelFont(42);
  h_residuals -> GetYaxis() -> SetTitleSize(0.17);
  h_residuals -> GetYaxis() -> SetTitleOffset(0.50);
  h_residuals -> GetYaxis() -> SetTitle("residuals");
  h_residuals -> GetYaxis() -> SetNdivisions(204);
  h_residuals -> GetXaxis() -> SetTitle("m_{ee} (GeV)");
  
  h_residuals -> SetMinimum(-5.);
  h_residuals -> SetMaximum(+5.);
  
  h_residuals -> Draw("hist");
  
  TGraphErrors* g_2sigma = new TGraphErrors();
  g_2sigma -> SetPoint(0,meeMin,0.);
  g_2sigma -> SetPointError(0,0.,2.);
  g_2sigma -> SetPoint(1,meeMax,0.);
  g_2sigma -> SetPointError(1,0.,2.);
  g_2sigma -> SetFillColor(kYellow);
  g_2sigma -> SetFillStyle(3001);
  g_2sigma -> Draw("E3,same");
  
  TGraphErrors* g_1sigma = new TGraphErrors();
  g_1sigma -> SetPoint(0,meeMin,0.);
  g_1sigma -> SetPointError(0,0.,1.);
  g_1sigma -> SetPoint(1,meeMax,0.);
  g_1sigma -> SetPointError(1,0.,1.);
  g_1sigma -> SetFillColor(kGreen);
  g_1sigma -> SetFillStyle(3001);
  g_1sigma -> Draw("E3,same");
  
  TF1* f_zero = new TF1("f_zero","0.",meeMin,meeMax);
  f_zero -> SetLineColor(kRed);
  f_zero -> SetLineWidth(2);
  f_zero -> Draw("same");
  
  h_residuals -> Draw("hist,same");
  
  
  p4 -> cd();
  p4 -> SetGridx();
  
  h_residualDistr -> GetXaxis() -> SetLabelSize(0.13);
  h_residualDistr -> GetXaxis() -> SetLabelFont(42);
  h_residualDistr -> GetXaxis() -> SetTitleSize(0.17);
  h_residualDistr -> GetXaxis() -> SetTitleOffset(1.00);
  h_residualDistr -> GetXaxis() -> SetTitle("");
  h_residualDistr -> GetYaxis() -> SetLabelSize(0.);
  h_residualDistr -> GetYaxis() -> SetTitleSize(0.17);
  h_residualDistr -> GetYaxis() -> SetTitleOffset(0.50);
  h_residualDistr -> GetYaxis() -> SetTitle("");
  h_residualDistr -> GetXaxis() -> SetNdivisions(204);
  h_residualDistr -> SetFillColor(kBlack);
  h_residualDistr -> SetFillStyle(3001);
  
  TF1* f_gaus = new TF1("f_gaus","[0]*exp(-1.*(x-[1])*(x-[1])/(2.*[2]*[2]))",-5.,5.);
  f_gaus -> SetParLimits(2,0.,10.);
  f_gaus -> SetParameters(1,0.,1.);
  h_residualDistr -> Fit("f_gaus","QRLS+");
  
  h_residualDistr -> Draw("");
  
  
  p3 -> cd();
  
  TLatex* latex;
  
  std::stringstream ss(label);
  std::string temp;
  int labelIt = 0;
  while(ss >> temp)
  {
    latex = new TLatex(0.10,0.85-0.10*labelIt,temp.c_str());
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.13);
    latex -> Draw("same");
    ++labelIt;
  }

  latex = new TLatex(0.10,0.20,Form("#chi^{2}/ndf = %1.2f",redChi2));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.13);
  latex -> Draw("");
  
  latex = new TLatex(0.10,0.10,Form("#mu = %1.2f #pm %1.2f",f_gaus->GetParameter(1),f_gaus->GetParError(1)));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.13);
  latex -> SetTextColor(kRed);
  latex -> Draw("same");
  
  latex = new TLatex(0.10,0.05,Form("#sigma = %1.2f #pm %1.2f",f_gaus->GetParameter(2),f_gaus->GetParError(2)));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.13);
  latex -> SetTextColor(kRed);
  latex -> Draw("same");
  
  gPad -> Update();  
  
  
  c -> Print((outFileName+"."+extension).c_str(),extension.c_str());
  delete c;
  delete h_residualDistr;
  
  outFile -> cd();
  
  h_mZ_MC[category] -> Write();
  h_mZ_DA[category] -> Write();
  h_mZFine_MC[category] -> Write();
  h_mZFine_DA[category] -> Write();
  h_mZRes_MC[category] -> Write();
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
