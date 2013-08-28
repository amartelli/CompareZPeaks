#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "ScaleEstimators.h"
#include "ZPeakFitUtils.h"

#include "TH1F.h"
#include "TRandom3.h"

#include <iostream>
#include <iomanip>



double meeFineMin = 60.;
double meeFineMax = 120.;
double precision = 0.01;



int main(int argc, char** argv)
{
  //---------------------------------------
  // Breit-Wigner and resolution parameters
  
  const float mZ = 91.1876;
  const float gammaZ = 2.4952;
  
  TF1* f_massReso_EB = new TF1("f_massReso_EB",crystalBallLowHigh,0.,2.,8);
  f_massReso_EB -> SetParameters(1.,1.00001e+00,1.17778e-02,1.40735e+00,1.58546e+00,7.43306e-01,2.13160e+01,1.00000e+00);
  
  TF1* f_massReso_EE = new TF1("f_massReso_EE",crystalBallLowHigh,0.,2.,8);
  f_massReso_EE -> SetParameters(1.,9.91089e-01,2.15497e-02,1.74378e+00,1.35601e+00,1.20988e+00,1.46226e+02,1.00000e+00);
  
  
  //-----------------
  // Other parameters
  
  double fraction = 0.6827;
  double mean,meanErr,min,max;
  
  int nTrials = 1000;
  int nEntries = 1000;
  if( argc > 1 ) nEntries = atoi(argv[1]);
  
  
  //------------
  // Begin loops
  
  TH1F* h_effectiveSigma_EB = new TH1F(Form("h_effectiveSigma_EB_%dEntries",nEntries),"",1000,0.,5.);
  TH1F* h_effectiveSigma_EE = new TH1F(Form("h_effectiveSigma_EE_%dEntries",nEntries),"",1000,0.,5.);
  
  for(int iTrial = 0; iTrial < nTrials; ++iTrial)
  {
    if( iTrial%50 == 0 )
      std::cout << ">>> trial " << iTrial << " / " << nTrials << " with " << nEntries << " entries each" << std::endl;
    
    TH1F* h_mZFine_EB = new TH1F("h_mZFine_EB","",int((meeFineMax-meeFineMin)/precision),meeFineMin,meeFineMax);
    TH1F* h_mZFine_EE = new TH1F("h_mZFine_EE","",int((meeFineMax-meeFineMin)/precision),meeFineMin,meeFineMax);
    
    for(int iEntry = 0; iEntry < nEntries; ++iEntry)
    {
      double m = gRandom -> BreitWigner(mZ,gammaZ);
      double r_EB = f_massReso_EB -> GetRandom();
      double r_EE = f_massReso_EE -> GetRandom();
      h_mZFine_EB -> Fill( m*r_EB );
      h_mZFine_EE -> Fill( m*r_EE );
    }
    
    FindSmallestInterval(mean,meanErr,min,max,h_mZFine_EB,fraction);
    h_effectiveSigma_EB -> Fill( 0.5*(max-min) );
    FindSmallestInterval(mean,meanErr,min,max,h_mZFine_EE,fraction);
    h_effectiveSigma_EE -> Fill( 0.5*(max-min) );
    
    delete h_mZFine_EB;
    delete h_mZFine_EE;
  }
  
  TFile* outFile = new TFile(Form("computeEffectiveSigmaErrors_%dEntries.root",nEntries),"RECREATE");
  h_effectiveSigma_EB -> Write();
  h_effectiveSigma_EE -> Write();
  outFile -> Close();
}
