#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "ScaleEstimators.h"

#include "TGraphErrors.h"
#include "TF1.h"

#include <iostream>
#include <iomanip>



double BW(double* x, double* par);
double GetBWIntegral(const double& fraction, TF1* f_BW, const double& prec = 0.001);



int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>> computeEffectiveSigma::usage: " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  
  //----------------------
  // Parse the config file
  
  parseConfigFile(argv[1]);
  
  std::string inputFile = gConfigParser -> readStringOption("Input::inputFile");
  
  float fraction = gConfigParser -> readFloatOption("Options::fraction");
  bool diagonalCatOnly = gConfigParser -> readBoolOption("Options::diagonalCatOnly");
  
  std::string outFilePath = gConfigParser -> readStringOption("Output::outFilePath");
  std::string outFileLabel = gConfigParser -> readStringOption("Output::outFileLabel");
  
  
  
  //---------------------------
  // Open input and output file
  std::cout << std::endl;
  std::cout << ">>> Open input and output file" << std::endl;
  
  TFile* inFile = TFile::Open(inputFile.c_str(),"READ");
  
  std::string outFileName = outFilePath;
  outFileName += "/computeEffectiveSigma";
  outFileName += "_" + outFileLabel;
  outFileName += Form("_%1.4f",fraction);
  
  TFile* outFile = TFile::Open((outFileName+".root").c_str(),"RECREATE");
  
  
  
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
  
  
  
  //------------------------
  // Compute effective sigma
  std::cout << std::endl;
  std::cout << ">>> Compute effective sigma" << std::endl;
  
  
  const float mZ = 91.188;
  const float gammaZ = 2.5;
  float xMin = -1710.;
  float xMax = +1890;
  TF1* f_BW = new TF1("f_BW",BW,xMin,xMax,2);
  f_BW -> SetParameters(mZ,gammaZ);  
  double seff_BW = GetBWIntegral(fraction,f_BW,0.0001);
  
  TF1* f_corr = new TF1("f_corr","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",0.,100.);
  f_corr -> SetParameters(1.30447e+01,-8.97976e+01,1.88517e+02,-1.58552e+02,4.77785e+01);
  
  TF1* f_err = new TF1("f_err","1.42734e/sqrt(x)",0.,1000000000);
  
  
  
  for(int R9Bin = 0; R9Bin < nR9Bins; ++R9Bin)
  {
    TGraphErrors* g_MC = new TGraphErrors();
    TGraphErrors* g_DA = new TGraphErrors();
    
    for(int etaBin = 0; etaBin < nEtaBins; ++etaBin)
    {
      std::string singleCatLabel = Form("eta%1.1f-%1.1f_R9%1.2f-%1.2f",etaBins.at(etaBin),etaBins.at(etaBin+1),R9Bins.at(R9Bin),R9Bins.at(R9Bin+1));
      std::string doubleCatLabel = singleCatLabel + "__" + singleCatLabel;
      std::cout << ">>>>>> processing category " << doubleCatLabel << std::endl;
      
      double mean_MC,meanErr_MC,min_MC,max_MC;
      std::string histoName_MC = Form("h_mZFine_MC_%s",(doubleCatLabel).c_str());
      TH1F* h_mZFine_MC = (TH1F*)( inFile->Get(histoName_MC.c_str()) );
      FindSmallestInterval(mean_MC,meanErr_MC,min_MC,max_MC,h_mZFine_MC,fraction);
      std::cout << ">>>>>>>>> MC: ["
                << std::fixed << std::setprecision(3) << std::setw(6) << min_MC
                << ","
                << std::fixed << std::setprecision(3) << std::setw(6) << max_MC
                << "]"
                << "   sigma_eff: "
                << std::fixed << std::setprecision(3) << std::setw(6) << 0.5*(max_MC-min_MC)
                << " ("
                << std::fixed << std::setprecision(3) << std::setw(5) << 0.5*(max_MC-min_MC)/mean_MC*100
                << "%)"
                << "   mean: "
                << std::fixed << std::setprecision(3) << std::setw(6) << mean_MC
                << " +/-"
                << std::fixed << std::setprecision(3) << std::setw(6) << meanErr_MC
                << std::endl;
      
      double mean_DA,meanErr_DA,min_DA,max_DA;
      std::string histoName_DA = Form("h_mZFine_DA_%s",(doubleCatLabel).c_str());
      TH1F* h_mZFine_DA = (TH1F*)( inFile->Get(histoName_DA.c_str()) );
      FindSmallestInterval(mean_DA,meanErr_DA,min_DA,max_DA,h_mZFine_DA,fraction);
      std::cout << ">>>>>>>>> DA: ["
                << std::fixed << std::setprecision(3) << std::setw(6) << min_DA
                << ","
                << std::fixed << std::setprecision(3) << std::setw(6) << max_DA
                << "]"
                << "   sigma_eff: "
                << std::fixed << std::setprecision(3) << std::setw(6) << 0.5*(max_DA-min_DA)
                << " ("
                << std::fixed << std::setprecision(3) << std::setw(5) << 0.5*(max_DA-min_DA)/mean_DA*100
                << "%)"
                << "   mean: "
                << std::fixed << std::setprecision(3) << std::setw(6) << mean_DA
                << " +/-"
                << std::fixed << std::setprecision(3) << std::setw(6) << meanErr_DA
                << std::endl;
      
      double x = 0.5 * (etaBins.at(etaBin) + etaBins.at(etaBin+1));
      double xErr = 0.5 * (etaBins.at(etaBin+1) - etaBins.at(etaBin));
      double seff_MC = 0.5*(max_MC-min_MC);
      double seff_DA = 0.5*(max_DA-min_DA);
      double BWInt_MC = f_BW->Integral(min_MC,max_MC)/fraction;
      double BWInt_DA = f_BW->Integral(min_DA,max_DA)/fraction;
      double yErr_MC = f_err->Eval(h_mZFine_MC->Integral()) * seff_MC / mean_MC;
      double yErr_DA = f_err->Eval(h_mZFine_DA->Integral()) * seff_DA / mean_DA;
      g_MC -> SetPoint(etaBin,x,sqrt(seff_MC*seff_MC-seff_BW*seff_BW)/mean_MC);
      g_MC -> SetPointError(etaBin,xErr,yErr_MC);
      g_DA -> SetPoint(etaBin,x,sqrt(seff_DA*seff_DA-seff_BW*seff_BW)/mean_DA);
      g_DA -> SetPointError(etaBin,xErr,yErr_DA);
    }
    
    outFile -> cd();
    g_MC -> Write(Form("g_MC_R9%1.2f-%1.2f",R9Bins.at(R9Bin),R9Bins.at(R9Bin+1)));
    g_DA -> Write(Form("g_DA_R9%1.2f-%1.2f",R9Bins.at(R9Bin),R9Bins.at(R9Bin+1)));
  }
  
}






double BW(double* x, double* par)
{
  return TMath::BreitWigner(x[0],par[0],par[1]);
}



double GetBWIntegral(const double& fraction, TF1* f_BW, const double& prec)
{
  double mZ = f_BW -> GetParameter(0);
  
  int step = 1;
  double integral = 0.;
  while( integral < fraction )
  {
    integral = f_BW -> Integral(mZ-step*prec,mZ+step*prec);
    ++step;
  }
  
  return step*prec;
}
