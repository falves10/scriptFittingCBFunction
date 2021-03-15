//#include "atlas/CommonHead.h"
#include <TTreeFormula.h>

#include "atlas/AtlasStyle.C"
#include "atlas/AtlasStyle.h"
#include "atlas/AtlasUtils.C"
#include "atlas/AtlasUtils.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TLorentzVector.h"

#include "TMath.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TSystem.h"

#include <math.h> // for “fabs”
#include <algorithm>
#include <map>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vector>
//RooFit include
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include <stdio.h>
#include <sys/stat.h>

#include "TMinuit.h"
using namespace RooFit ;
using namespace std;


class FitInfo
{ 
public:
    int eta; 
    int energy; 
    int bias;
    double constant; 
    double mean;
    double sigma;
    double alpha; 
    double n;
    double error;
    bool broken;
    FitInfo(int eta_, int energy_, int bias_, double constant_, double mean_, double sigma_, double alpha_, double n_, double error_):
    eta(eta_),
    energy(energy_),
    bias(bias_),
    constant(constant_),
    mean(mean_),
    sigma(sigma_),
    alpha(alpha_),
    n(n_),
    error(error_),
    broken(false)
    {}
};

//linear function used to fix broken data
class LinearFitFunc
{
public:
    LinearFitFunc(double alpha, double beta): _alpha(alpha), _beta(beta){}
    LinearFitFunc(const LinearFitFunc& a)
    {
        this->_alpha = a._alpha;
        this->_beta = a._beta;
    }
    
    double Func(double x)
    {
        return _beta * x + _alpha;
    }
private:
    double _alpha;
    double _beta;
};

std::vector<std::shared_ptr<FitInfo>> gInfo;
std::vector<std::shared_ptr<FitInfo>> gMCInfo;
std::vector<std::shared_ptr<FitInfo>> gDataInfo;
typedef std::pair<int,int> BinEnergyPair;
std::map<BinEnergyPair, std::vector<std::shared_ptr<LinearFitFunc>>> gFixFuncMap;
std::map<int, std::vector<std::shared_ptr<LinearFitFunc>>> gFixFuncMCMap;
std::map<int, std::vector<std::shared_ptr<LinearFitFunc>>> gFixFuncDataMap;

void Initialize();
void LoadFitInfo(string path);
vector<std::shared_ptr<FitInfo>> GetFitInfo(int iEtaBin, int iEnergyzBin, bool excludingBiasMinusOne);
vector<std::shared_ptr<FitInfo>> GetFitInfoMCSimulation(int iEtaBin);
vector<std::shared_ptr<FitInfo>> GetFitInfoExperimentalData(int iEtaBin);
vector<std::shared_ptr<LinearFitFunc>> GenerateFitFuncBiasMinusOne(const vector<std::shared_ptr<FitInfo>>& info);
vector<std::shared_ptr<LinearFitFunc>> GenerateFitFunc(const vector<std::shared_ptr<FitInfo>>& info);
void MarkBrokenFit(vector<std::shared_ptr<FitInfo>>& v,  double threshold);
void FindBrokenBias(vector<FitInfo>& v);


//int m_nEtaBin(13);
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.37, 1.52,
//     1.8, 1.9, 2.1, 2.3, 2.47});

int m_nEtaBin(16);
//int m_nEtaBin(14);
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.37, 1.52,1.8, 1.9, 2.1, 2.3, 2.47});
std::vector<double> m_etaBinning({0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30,1.35,1.40});
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.775, 0.825, 1., 1.2, 1.37, 1.52,1.8, 1.9, 2.1, 2.3, 2.47});
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.7, 0.825, 1., 1.2, 1.37, 1.52,1.8, 1.9, 2.1, 2.3, 2.47});


//int m_nEnergyBin(10);
int m_nEnergyBin(10);
//std::vector<double> m_energyBinning({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.0});
std::vector<double> m_energyBinning({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.0});


bool calibrated(true);//Adele true
float lowcut = 87.;
float highcut = 94.5;

bool doMaterial = false;

void PrintProgressBar( int index, int total )
{
  if( index%200000 == 0 )
    {
      TString print_bar = " [";
      for( int bar = 0; bar < 20; bar++ )
	{
	  double current_fraction = double(bar) / 20.0;
	  if( double(index)/double(total) > current_fraction )
	    print_bar.Append("/");
	  else
	    print_bar.Append(".");
	}
      print_bar.Append("] ");
      std::cout << print_bar << 100.*(double(index)/double(total)) << "%\r" << std::flush;
    }
}



void infile2chain(TString _infilelist, TChain *&_fchain, TString chainname)
{
  _fchain = new TChain(chainname);
  ifstream infile(_infilelist, ios::in);
  string line;
  while (getline(infile, line)){
    _fchain->Add(line.c_str());
  }
  infile.close();
}

//convolution function for the mass
std::vector<float> convolution(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, double input_width,double input_sigma, TString path, int m_layer);

std::vector<float> convolution_eOverp(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, TString path, int m_layer, TString campaign, std::ofstream& file_brokenfit);

//methods to save plots
void mySave(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, int iEnergyBin,TString vartagname, TString varname,TString path, int m_layer);
void mySave2D(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, TString vartagname, TString varname, TString path,TString parameter, int m_layer);

vector<double> fixingOnBadFits( int iEtaBin, int iEnergyBin, int iBias );
bool checkedFailedFit( int iEtaBin, int iEnergyBin, int iBias );
void simpleOnBadFits( int iEtaBin, int iEnergyBin, int iBias );
vector<double> testeFunctional(int iEtaBin, int iEnergyBin, int iBias);
bool throwKinematics(int iEtaBin, int iEnergyBin, int iBias);

/// Get the bin number of a given eta value
unsigned int binNbOfEta(double eta) 
{
  double localEta = eta;
  localEta = fabs(eta);
  unsigned int bin = 0;
  while (bin < m_nEtaBin && localEta > m_etaBinning[bin+1])
    bin++;
  return bin;
} // end binNbOfEta


/// Get the bin number of a given E1/E2 value
unsigned int binNbOfEnergy(double E1E2) 
{
  
  double localEnergy = E1E2;
  localEnergy = fabs(E1E2);
  unsigned int bin = 0;
  while (bin < m_nEnergyBin && localEnergy > m_energyBinning[bin+1])
    bin++;
  return bin;
} // end binNbOfE1E2




//MAIN
//int invariantMassFitCal2(TString input) {
int invariantMassFitOnly(TString campaign,TString layer, TString variable) {
 
  SetAtlasStyle();  
   
  vector<TString> varname;
  varname.push_back(variable);
  
  vector<TString> vartagname;
  if(variable.Contains("mass"))
    vartagname.push_back("m_{ee} [GeV]");
  else if (variable.Contains("eOverp"))
    vartagname.push_back("E/p");

  int m_layer = atoi(layer);
  std::cout<<"m_layer = "<<m_layer<<std::endl;
    
  for(int ivar = 0; ivar < varname.size(); ivar++){

    std::cout<<"variable = "<<varname[ivar]<<std::endl;

    //TString path = Form("/publicfs/atlas/atlasnew/higgs/hgg/donofrioadele/EG_condor/EGAMMAStep2Code/output_"+varname[ivar]+"_"+campaign+"_layer%i_%f_%f/",m_layer,lowcut,highcut);
      
    TString path = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/output_twoLeptons_Nominal/";
    gSystem->mkdir(path, kTRUE);

    std::cout<<"path_out = "<<path<<std::endl;

    ofstream myfile_brokenfits;
    myfile_brokenfits.open ("out_fit.dat");  
      
      
    TH1F* map_histog[17];
    TH1F* vhisto2D_data[17];
    TH1F* vhisto2DMean_data[17];
    TH1F* vhisto2DSigma_data[17];
    
    std::map<pair<int, int>, map<int,TH1F*>> m_histo; 
    
    int biasSize = 16;
    //int biasSize = 2;
    
    int binning; float lowedge; float highedge; float lowcut; float highcut;// mass
    if(varname[ivar].Contains("mass")) {binning = 50; lowedge = 80; highedge = 100; lowcut = 80; highcut = 100;}// mass: 80-100
    if(varname[ivar].Contains("m_eOverp")) {binning = 50; lowedge = 0.6; highedge = 2; lowcut = 0.6; highcut = 2;}// E/p
    
    unsigned int maxEta = 0;
    
    maxEta = m_nEtaBin;
    //get input files
    //TString path_in = Form("/publicfs/atlas/atlasnew/higgs/hgg/donofrioadele/EG_condor/EGAMMAStep1Code/adele_ntuple/calibrated_"+varname[ivar]+"/Histograms/layer%i/",m_layer);
      
    TString path_in = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/";

    
     
      
      
    float array[11] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.0};
    
    TH1F* histo2D_mc = new TH1F ("mass_etabin_mc", "", m_nEnergyBin,array);
    TH1F* histo2D_data = new TH1F ("mass_etabin_data", "", m_nEnergyBin,array);
    //Mean mass method
    TH1F* histo2DMean_mc = new TH1F ("mass_etabin_mcMean", "", m_nEnergyBin,array);
    TH1F* histo2DMean_data = new TH1F ("mass_etabin_dataMean", "", m_nEnergyBin,array);
    //CB Parameters values
    TH1F* histo2DSigma_data = new TH1F ("mass_etabin_dataSigma", "", m_nEnergyBin,array);
    TH1F* histo2DSigma_mc = new TH1F ("mass_etabin_mcSigma", "", m_nEnergyBin,array);    
    
    histo2D_mc->Reset();
    histo2D_data->Reset();
    histo2DMean_mc->Reset();
    histo2DMean_data->Reset();
    histo2DSigma_mc->Reset();
    histo2DSigma_data->Reset();
    
    std::cout<<"biasSize = "<<biasSize<<std::endl;

    //weight_histo_eta7_energy9_bias3_m_eOverp_mergedOutput.root   
      
    TString forMaterial = "_FMX";
    TFile *f_in;  
      
    if(doMaterial)
    {
        f_in = new TFile(path_in+"weight_histo_"+campaign+"_"+varname[ivar]+forMaterial+".root");
        if ( !(f_in->IsOpen()) )continue; 
        
        std::cout<<"path_in = "<<path_in+"weight_histo_"+campaign+"_"+varname[ivar]+forMaterial+".root"<<std::endl;
    }
    else
    {
        f_in = new TFile(path_in+"weight_histo_"+campaign+"_"+varname[ivar]+".root");
        if ( !(f_in->IsOpen()) )continue; 
        
        std::cout<<"path_in = "<<path_in+"weight_histo_"+campaign+"_"+varname[ivar]+".root"<<std::endl;
    }
      

    

    for (unsigned int iEtaBin = 0; iEtaBin < maxEta; iEtaBin++) {
        
      //if( iEtaBin != 0 ) continue;
      cout<<" Fabio_lep_eta " << iEtaBin << endl;
      for(int i = 2; i < biasSize+2 ; i++){//Adele
          
	   vhisto2D_data[i-2] = new TH1F (Form("mass_etabin%d_data_bias%i",iEtaBin,i), "", m_nEnergyBin,array);
	   vhisto2DMean_data[i-2] = new TH1F (Form("massMean_etabin%d_data_bias%i",iEtaBin,i), "", m_nEnergyBin,array);
	   vhisto2DSigma_data[i-2] = new TH1F (Form("mass_etabin%d_data_biasSigma%i",iEtaBin,i), "", m_nEnergyBin,array);
          
          
      }
      
      /*for(int i = 0; i < 17; i++)
      {
          cout << vhisto2D_data[i]->GetName() << endl;
      }*/  
        
      for (unsigned int iEnergyBin = 0; iEnergyBin < m_nEnergyBin; iEnergyBin++) {
          
        //if(iEnergyBin != 1) continue;
          
	map_histog[1] = new TH1F(Form("data_%d_%d",iEtaBin,iEnergyBin), "",    binning,  lowedge, highedge);
	map_histog[0] = new TH1F(Form("mc_%d_%d",iEtaBin,iEnergyBin), "", binning,  lowedge, highedge);
	
	map_histog[0] = (TH1F*) f_in->Get(Form("mc_eta_%d_energy_%d",iEtaBin,iEnergyBin));
	map_histog[1] = (TH1F*) f_in->Get(Form("data_eta_%d_energy_%d",iEtaBin,iEnergyBin));
       
	if(map_histog[0]->Integral()==0) continue;
	if(map_histog[1]->Integral()==0) continue;
	
	for(int i = 0; i < biasSize ; i++)
    {//Adele
	    map_histog[i+2] = (TH1F*) f_in->Get(Form("data_eta_%d_energy_%d_bias%d",iEtaBin,iEnergyBin,i));
	    map_histog[i+2]->Sumw2();
	  }
	  
	  histo2D_mc->Sumw2();
	  histo2D_data->Sumw2();
	  
	  ///////////////////////////////////////////////////////////////////////////////
	  float mean_mass_mc = -1.;
	  float mean_mass_data = -1.;
	  float emean_mass_mc = -1.;
	  float emean_mass_data = -1.;
	  float width_mass_mc = 2.49;
	  float width_mass_data = -1.;
	  float ewidth_mass_mc = -1.;
	  float ewidth_mass_data = -1.;
	  float sigma_mass_mc = 2.6;
	  float sigma_mass_data = -1.;
	  float esigma_mass_mc = -1.;
	  float esigma_mass_data = -1.;
	  
	  std::vector<float> result_data;
	  std::vector<float> result_mc;

	  if(varname[ivar].Contains("mass"))
	    result_mc = convolution(map_histog[0],-1,iEtaBin,iEnergyBin,1,map_histog[0]->GetMean(),width_mass_mc,sigma_mass_mc,path, m_layer);
	  else if(varname[ivar].Contains("eOverp"))
      {
          if(map_histog[0]->GetEntries() < 1000) break; //for new binning, data/MC
          //if(map_histog[0]->GetEntries() <= 00) break; 
          else result_mc = convolution_eOverp(map_histog[0],-1,iEtaBin,iEnergyBin,1,map_histog[0]->GetMean(),path, m_layer, campaign, myfile_brokenfits);
      }
           
	
	  mean_mass_mc = result_mc.at(0);
	  emean_mass_mc = result_mc.at(1);
          
	  if(varname[ivar].Contains("mass"))
      {
	    sigma_mass_mc = result_mc.at(3);
	    esigma_mass_mc = result_mc.at(4);
	  }
          
	  if(varname[ivar].Contains("mass"))
	    result_data = convolution(map_histog[1],-1,iEtaBin,iEnergyBin,0,mean_mass_mc,width_mass_mc,sigma_mass_mc,path, m_layer);
	  else if(varname[ivar].Contains("eOverp"))
      {
          if(map_histog[1]->GetEntries() < 1000 ) break;////for new binning, data/MC
          //if(map_histog[1]->GetEntries() <= 100) break;
          else result_data = convolution_eOverp(map_histog[1],-1,iEtaBin,iEnergyBin,0,mean_mass_mc,path, m_layer, campaign, myfile_brokenfits);
      }
	      
	  
          
      mean_mass_data = result_data.at(0);
	  emean_mass_data = result_data.at(1);
	  
      if(varname[ivar].Contains("mass"))
      {
	    sigma_mass_data = result_data.at(3);
	    esigma_mass_data = result_data.at(4);
	  }
          
	  histo2D_mc->SetBinContent(iEnergyBin+1,mean_mass_mc);
	  histo2D_data->SetBinContent(iEnergyBin+1,mean_mass_data);
	  histo2D_mc->SetBinError(iEnergyBin+1,emean_mass_mc);
	  histo2D_data->SetBinError(iEnergyBin+1,emean_mass_data);

	  if(varname[ivar].Contains("mass")){
	    //Mean mass method
	    histo2DMean_mc->SetBinContent(iEnergyBin+1,map_histog[0]->GetMean());
	    histo2DMean_mc->SetBinError(iEnergyBin+1,map_histog[0]->GetMeanError());
	    histo2DMean_data->SetBinContent(iEnergyBin+1,map_histog[1]->GetMean());
	    histo2DMean_data->SetBinError(iEnergyBin+1,map_histog[1]->GetMeanError());
	    //Plot all the CB+BW function parameters
	    histo2DSigma_data->SetBinContent(iEnergyBin+1,result_data.at(3));
	    histo2DSigma_data->SetBinError(iEnergyBin+1,result_data.at(4));
	    histo2DSigma_mc->SetBinContent(iEnergyBin+1,result_mc.at(3));
	    histo2DSigma_mc->SetBinError(iEnergyBin+1,result_mc.at(4));
	  }

	  for(int i = 2; i < biasSize+2 ; i++)
      {//Adele
          
        
          
	    float vmean_mass_data = -1;
	    float vemean_mass_data = -1;
	    
	    std::vector<float> vresult;
	    vresult.clear();
	    int binmax = map_histog[i]->GetMaximumBin(); 
	    double x = map_histog[i]->GetXaxis()->GetBinCenter(binmax);
        //map_histog[i]->Print("all");
	    if(varname[ivar].Contains("mass"))
	      vresult = convolution(map_histog[i],i,iEtaBin,iEnergyBin,0,mean_mass_mc,width_mass_mc,sigma_mass_mc,path, m_layer);
          
	    if(varname[ivar].Contains("eOverp"))
            //if(iEtaBin == 0 && i == 2)
            {
                if(map_histog[i]->GetEntries() < 1000) break;//for new binning, data/MC
                //if(map_histog[i]->GetEntries() < 100) break;
                else vresult = convolution_eOverp(map_histog[i],i,iEtaBin,iEnergyBin,0,mean_mass_mc,path, m_layer, campaign, myfile_brokenfits);  
                
                vmean_mass_data = vresult.at(0);
	             vemean_mass_data = vresult.at(1);
            }
	        
	    //vmean_mass_data = vresult.at(0);
	    //vemean_mass_data = vresult.at(1);
        cout<< vmean_mass_data << " " << vemean_mass_data << endl;
        
	    vhisto2D_data[i-2]->SetBinContent(iEnergyBin+1,vresult.at(0));
	    vhisto2D_data[i-2]->SetBinError(iEnergyBin+1,vresult.at(1));
        vhisto2D_data[i-2]->Print("all");
	    if(varname[ivar].Contains("mass")){
	      //Mean mass method
	      vhisto2DMean_data[i-1]->SetBinContent(iEnergyBin+1,map_histog[i]->GetMean());
	      vhisto2DMean_data[i-1]->SetBinError(iEnergyBin+1,map_histog[i]->GetMeanError());
	      //Plot all the CB+BW function parameters
	      vhisto2DSigma_data[i-1]->SetBinContent(iEnergyBin+1,vresult.at(3));
	      vhisto2DSigma_data[i-1]->SetBinError(iEnergyBin+1,vresult.at(4));
	    }
	  }
          
	  mySave(map_histog[0],map_histog[1],-1,iEtaBin,iEnergyBin,vartagname[ivar],varname[ivar],path, m_layer);
	  for(int i = 0; i < biasSize ; i++)
      {//Adele
	    mySave(map_histog[0],map_histog[i+2],i+2,iEtaBin,iEnergyBin,vartagname[ivar],varname[ivar],path, m_layer);
	  }
          
	}//end loop energy bins
      
    
	mySave2D(histo2D_mc,histo2D_data,-1,iEtaBin,vartagname[ivar],varname[ivar],path,"", m_layer);
	if(varname[ivar].Contains("mass")){
	  mySave2D(histo2DMean_mc,histo2DMean_data,-1,iEtaBin,vartagname[ivar],varname[ivar],path,"Mean", m_layer);
	  mySave2D(histo2DSigma_mc,histo2DSigma_data,-1,iEtaBin,vartagname[ivar],varname[ivar],path,"Sigma", m_layer);
	}
    
	for(int i = 0; i < biasSize ; i++)
    {//Adele
      
	  mySave2D(histo2D_mc,vhisto2D_data[i],i+2,iEtaBin,vartagname[ivar],varname[ivar],path,"", m_layer);
	    if(varname[ivar].Contains("mass")){
	      mySave2D(histo2DMean_mc,vhisto2DMean_data[i],i+2,iEtaBin,vartagname[ivar],varname[ivar],path,"Mean", m_layer);
	      mySave2D(histo2DSigma_mc,vhisto2DSigma_data[i],i+2,iEtaBin,vartagname[ivar],varname[ivar],path,"Sigma", m_layer);
	    }
	}
	//saving 2D histos
      }//end loop eta bins
      
    delete histo2D_mc;
    delete histo2D_data;
    delete histo2DSigma_mc;
    delete histo2DSigma_data;
    delete histo2DMean_mc;
    delete histo2DSigma_mc;
      
    myfile_brokenfits.close();
      
  }//end loop variables
  return 0;
}



std::vector<float> convolution(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, double input_width, double input_sigma, TString path, int m_layer){	

  std::vector<float> parameters;
  parameters.clear();

  //change the range of the plots to get the mean!!!
  RooRealVar x( "x", "x", 87., 94.5 );//84,98//80-100
  x.setBins(10000,"cache") ;
  x.setMin("cache",64.) ;
  x.setMax("cache",118.) ;
  
  //std::cout<<"mean_bw = "<<mean_bw<<std::endl;
  // Breit-Wigner                        
  RooRealVar m0( "m0", "m0", mean_bw, 87., 94.5);//80-100
  RooRealVar width( "width", "width", 2.49);//2.49,1.,4.
  if(iperiod == 0){
      width.setConstant(kTRUE);
  }
  RooBreitWigner bw( "bw", "bw", x, m0, width);

  // Crystal-Ball                                                                                                                         
  RooRealVar mean( "mean", "mean", 0. );
  RooRealVar sigma( "sigma", "sigma", input_sigma, 1., 5.);//2.6,1.,5.
  if(iperiod == 0){
    sigma.setConstant(kTRUE);
  }
  RooRealVar alpha( "alpha", "alpha", 1.3 );
  RooRealVar n( "n", "n", 5.1 );
  RooCBShape cb( "cb", "cb", x, mean, sigma, alpha, n );
  // convolution                                                                                                                         
  RooFFTConvPdf pdf_sig( "pdf_sig", "pdf_sig", x, bw, cb );
  //add polynomial bkg
  // Build Chebychev polynomial p.d.f.
  RooRealVar a0("a0","a0",0.5,0.,1.) ;
  RooRealVar a1("a1","a1",0.,0.,1.) ;
  RooChebychev bkg1("bkg1","bkg1",x,RooArgSet(a0,a1));
  // Sum the composite signal and background
  RooRealVar bkgfrac("bkgfrac","fraction of background",0.05,0.,1.);
  RooAddPdf  pdf("pdf","pdf",RooArgList(bkg1,pdf_sig),bkgfrac) ;

  RooDataHist histo("histo","histo",x,Import(*map_histog));

  if(iperiod == 0)
    pdf.fitTo(histo,SumW2Error(kTRUE)) ;
  if(iperiod == 1)
    pdf_sig.fitTo(histo,SumW2Error(kTRUE)) ;

  TCanvas canv( "canv", "canv", 800., 600. );
  RooPlot* frame1 = x.frame(Bins(100),Title("Convolution of a Breit-Wigner and a Crystal-Ball, Chebychev pol. bkg")) ;
  histo.plotOn(frame1,Name("Data")) ;
  if(iperiod == 0){
    pdf.plotOn(frame1,Name("pdf"),LineColor(kRed)) ;
    pdf.paramOn(frame1,Layout(0.60));
    pdf.plotOn(frame1,Components("bkg1"),LineStyle(kDotted),LineColor(kBlue));
  }
  if(iperiod == 1){
    pdf_sig.plotOn(frame1,Name("pdf"),LineColor(kRed)) ;
    pdf_sig.paramOn(frame1,Layout(0.60));
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",800,600) ;
  canvas->cd() ; 
  TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
  TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
  pad1->Draw();pad2->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.16);
  pad2->SetBottomMargin(0.24);
  frame1->GetYaxis()->SetTitleOffset(1.4) ;
  frame1->GetXaxis()->SetTitle((Form("m_{ee} [GeV], #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i ", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));
  frame1-> Draw();
  pad1->Modified();
  pad1->RedrawAxis();
  pad1->Update();
  pad2->cd();

  TF1* tf1_model_data = pdf.asTF(x);
  TF1* tf1_model = pdf_sig.asTF(x);

  TH1* clone_data = (TH1*)histo.createHistogram("clone_data",x,Binning(50));
  RooDataHist *pdfHisto_data = pdf.generateBinned(x,1000000);
  TH1* clone_fit_data = (TH1*)pdfHisto_data->createHistogram("clone_fit_data",x,Binning(50));
  clone_fit_data->Scale(clone_data->Integral()/clone_fit_data->Integral());

  RooDataHist *pdfHisto = pdf_sig.generateBinned(x,1000000);
  TH1* clone_fit = (TH1*)pdfHisto->createHistogram("clone_fit",x,Binning(50));
  clone_fit->Scale(clone_data->Integral()/clone_fit->Integral());

  if(iperiod == 0)
    clone_data->Divide(clone_fit_data);

  else if(iperiod == 1)
    clone_data->Divide(clone_fit);

  clone_data->GetXaxis()->SetTitle((Form("m_{ee} [GeV], #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));
  clone_data->GetYaxis()->SetTitle("DATA / FIT");
  clone_data->GetXaxis()->SetRangeUser(lowcut, highcut);//80-100
  clone_data->GetYaxis()->SetRangeUser(0.7,1.3);
  clone_data->GetXaxis()->SetLabelSize(0.1);
  clone_data->GetYaxis()->SetLabelSize(0.08);
  clone_data->GetXaxis()->SetTitleSize(0.08);
  clone_data->GetYaxis()->SetTitleSize(0.09);
  clone_data->GetYaxis()->SetTitleOffset(0.6);
  clone_data->GetXaxis()->SetTitleOffset(1.2);
  clone_data->Draw("E1");
  pad2->Modified();
  pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();

  TString output_folder = path+Form("FitPlots/layer%i/",m_layer);
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } 
  else {
    gSystem->mkdir(path+Form("FitPlots/layer%i/",m_layer), kTRUE);
  }

  canvas->SaveAs(path+Form("FitPlots/layer%i/weight_period%i_invariantMass_eta%i_energy%i_bias%i.pdf",m_layer,iperiod,iEtaBin,iEnergyBin,iBias));
  canvas->SaveAs(path+Form("FitPlots/layer%i/weight_period%i_invariantMass_eta%i_energy%i_bias%i.root",m_layer,iperiod,iEtaBin,iEnergyBin,iBias));

   
  double modelMean_data = tf1_model_data->GetMaximumX();
  double modelMean = tf1_model->GetMaximumX();

  if(iperiod == 0)
    parameters.push_back(modelMean_data);//0
  if(iperiod == 1)
    parameters.push_back(modelMean);//0
  parameters.push_back(m0.getError());//1
  parameters.push_back(m0.getVal());//2
  parameters.push_back(sigma.getVal());//3
  parameters.push_back(sigma.getError());//4
  parameters.push_back(mean.getVal());//5
  parameters.push_back(mean.getError());//6
  parameters.push_back(alpha.getVal());//7
  parameters.push_back(alpha.getError());//8
  parameters.push_back(n.getVal());//9
  parameters.push_back(n.getError());//10
  parameters.push_back(width.getVal());//11
  parameters.push_back(width.getError());//12

  cout << "width = "<<width.getVal() << endl;

  delete frame1;
  delete clone_data;
  delete clone_fit;
  delete pad1;
  delete pad2;
  delete canvas;

  return parameters;

}


//////////////////////////////////////////eOverP////////////////////////////////////////
 std::vector<float> convolution_eOverp(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, TString path, int m_layer, TString campaign, std::ofstream& file_brokenfit)
{
     cout<<" ETA " << iEtaBin << " ENERGY " << iEnergyBin << " BIAS "<< iBias << endl;
     cout<<" ENTRIES " << map_histog->GetEntries() << endl;
   

    std::vector<TString> parname_ws    = {"Constant", "Mean","Sigma","Alpha", "N"};
    std::vector<float> parameters;
    parameters.clear();
    
TF1 *f2 = new TF1("f2","crystalball",0.9,1.3); //Nominal fit range 
//TF1 *f2 = new TF1("f2","crystalball",0.95,1.2); //Alternative fit range 
     
f2->SetLineColor(2);
bool teste = false;
bool noParms = false; //for other fitrange
bool isNomFitRange = false;
bool isPP0 = false;
bool isNomElectronID = true; //used false for etabin4 data 2017 agains pp0 config
     
bool globalRun = false;
  

if(globalRun)
{
    if( iperiod == 0 )
    {
        if(campaign == "mc16")
        {
            if(iEtaBin == 0 || iEtaBin == 1 || iEtaBin == 2 || iEtaBin == 3)
            {
                /*
                        if( iEnergyBin == 6)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else if( iEnergyBin == 7 || iEnergyBin == 8 || iEnergyBin == 9 )
                        {
                            teste = true;
                            f2->SetParameters( 6.34042e+00 ,  1.04745e+00  , 5.78203e-02 , -5.96818e-01 , 1.66073e+00);
                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00); 
                */ 
                /*if( (iEnergyBin == 2 && iBias == -1)  || (iEnergyBin == 3 && iBias == -1) || (iEnergyBin == 4 && iBias == -1)  )
                {
                    teste = true;
                    f2->SetParameters( 1.11068e+01  ,  9.99187e-01  , 5.03072e-02   ,  -5.77218e-01 , 4.34526e+00);
                }*/
                if( iEnergyBin == 5  )
                {
                    teste = true;
                    f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                }
                else if( iEnergyBin == 6 || iEnergyBin == 7 )
                {
                    teste = true;
                    f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                }
                else if( iEnergyBin == 8 || iEnergyBin == 9 )
                {
                    teste = true;
                    f2->SetParameters( 8.23483e+00  , 9.89187e-01  , 3.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                }
                else
                {
                    teste = true;
                    f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                }


            }
            else if(iEtaBin == 3 )
            {
                    /*
                    if( iEnergyBin == 4)
                        {
                            if(iBias == 9 || iBias == 11 || iBias == 15 )
                            {
                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                            }else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                        }
                        else if( iEnergyBin  == 5 || iEnergyBin  == 6 )
                        {

                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );

                        }
                        else if ( iEnergyBin  == 7 || iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else if ( iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        */

                /*if( ( iEnergyBin == 3 && iBias == -1 ) || ( iEnergyBin == 4 && iBias == -1 )  )
                {
                    teste = true;
                    f2->SetParameters( 1.11068e+01  ,  9.99187e-01  , 5.03072e-02   ,  -5.77218e-01 , 4.34526e+00);
                }*/
                if ( iEnergyBin  == 9)
                {
                    teste = true;
                    f2->SetParameters( 6.32709e-02  , 1.02743e+00  , 6.28142e-02   ,  -6.33010e-01  ,  1.03690e+00 );
                }
                else
                {
                    teste = true;
                    f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
            }
            else if(iEtaBin == 4 || iEtaBin == 5 )
            {
                    /*
                    if( iEnergyBin == 4)
                        {
                            if(iBias == 9 || iBias == 11 || iBias == 15 )
                            {
                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                            }else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                        }
                        else if( iEnergyBin  == 5 || iEnergyBin  == 6 )
                        {

                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );

                        }
                        else if ( iEnergyBin  == 7 || iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else if ( iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        */


                        /*
                        if( ( iEnergyBin == 2 && iBias == -1 ) || ( iEnergyBin == 3 && iBias == -1 ) || ( iEnergyBin == 4 && iBias == -1 ) || ( iEnergyBin == 5 && iBias == -1 ) || ( iEnergyBin == 6 && iBias == -1 )   )
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  9.99187e-01  , 5.03072e-02   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        */
                        if ( iEnergyBin  == 7  && iBias == 4)
                        {
                            teste = true;
                            f2->SetParameters( 6.32709e-02  , 1.02743e+00  , 6.28142e-02   ,  -6.33010e-01  ,  1.03690e+00 );
                        }
                        else if ( iEnergyBin  == 8 && iBias == 8)
                        {
                            teste = true;
                            f2->SetParameters( 6.32709e-02  , 1.02743e+00  , 1.28142e-01   ,  -6.33010e-01  ,  1.03690e+00 );
                        }
                        else if ( iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters( 6.32709e-02  , 1.02743e+00  , 6.28142e-02   ,  -6.33010e-01  ,  1.03690e+00 );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 6 )
            {
                    /*
                    if( iEnergyBin == 4)
                        {
                            if(iBias == 9 || iBias == 11 || iBias == 15 )
                            {
                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                            }else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                        }
                        else if( iEnergyBin  == 5 || iEnergyBin  == 6 )
                        {

                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );

                        }
                        else if ( iEnergyBin  == 7 || iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else if ( iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        */

                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters(4.05517e-01 , 1.01206e+00 , 5.83740e-02  , -4.39952e-01 , 2.26315e+00  );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 7 )
            {


                        if ( iEnergyBin  == 0)
                        {


                            {
                                teste = true;
                                f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                            }

                        }
                        else if ( iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters(1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00);
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 8 )
            {


                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 4 || iEnergyBin  == 5 || iEnergyBin  == 6)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00 , 1.00662e+00 , 4.42452e-02  , -5.86730e-01  , 1.10299e+00 );
                        }
                        else if ( iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00, 1.10356e+00 , 1.00000e-02  , -1.43864e+00   ,1.00000e+00  );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 9 )
            {


                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 4 || iEnergyBin  == 5 || iEnergyBin  == 6)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00 , 1.00662e+00 , 4.42452e-02  , -5.86730e-01  , 1.10299e+00 );
                        }
                        else if ( iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00, 1.10356e+00 , 1.00000e-02  , -1.43864e+00   ,1.00000e+00  );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 10 )
            {


                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 4 || iEnergyBin  == 5 || iEnergyBin  == 6)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 6.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 7 ||  iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 6.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 11 )
            {


                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 4 || iEnergyBin  == 5 || iEnergyBin  == 6)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 6.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 7)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 1.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00 , 9.92017e-01, 5.48004e-02   , -1.39785e+00  , 1.00000e+00 );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 12 )
            {


                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 4 || iEnergyBin  == 5 || iEnergyBin  == 6)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 6.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 7)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 1.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00 , 9.92017e-01, 5.48004e-02   , -1.39785e+00  , 1.00000e+00 );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if( iEtaBin == 13 )
                    {


                        if ( iEnergyBin  == 0)
                        {
                            teste = true;
                            f2->SetParameters( 1.64301e+00  , 9.91284e-01 , 4.31344e-02   , -1.02052e+00 , 2.58542e+00 );
                        }
                        else if ( iEnergyBin  == 4 || iEnergyBin  == 5 || iEnergyBin  == 6)
                        {
                            teste = true;
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 6.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 7)
                        {
                            teste = true;
                            //f2->SetParameters(7.33553e+00 , 1.03037e+00, 1.43383e-02   , -6.84461e-01 , 2.03392e+00);
                            f2->SetParameters(7.33553e+00 , 1.03037e+00, 6.43383e-02   , -6.84461e-01 , 2.03392e+00);
                        }
                        else if ( iEnergyBin  == 8)
                        {
                            teste = true;
                            f2->SetParameters(1.00000e+00 , 9.92017e-01, 5.48004e-02   , -1.39785e+00  , 1.00000e+00 );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                        }

                        //else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            /*
            else if (iEtaBin == 1)
            {
                        if( iEnergyBin == 7 || iEnergyBin == 8 || iEnergyBin == 9    )
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00); 
            }*/
            /*
            else if (iEtaBin == 2)
            {
                        if( iEnergyBin == 7 || iEnergyBin == 8 || iEnergyBin == 9 )
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00); 
                    }
            else if(iEtaBin == 3)
            {
                        if( iEnergyBin == 0 && iBias == 3)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else if( iEnergyBin == 0 && iBias == 5)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else if(iEnergyBin == 8 && iBias == 2)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else if(iEnergyBin == 8 && iBias == 4)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else if(iEnergyBin == 9 && iBias == 2)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else if(iEnergyBin == 9 && iBias == 4)
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if(iEtaBin == 4)
            {
                        if( iEnergyBin == 4)
                        {
                            if(iBias == 9 || iBias == 11 || iBias == 15 )
                            {
                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );
                            }else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                        }

                        else if( iEnergyBin  == 5 || iEnergyBin  == 6 )
                        {

                                teste = true;
                                f2->SetParameters( 1.53140e+01  , 1.03003e+00  , 6.48004e-02 , -5.79820e-01 , 1.85537e+00 );

                        }
                        else if ( iEnergyBin  == 7 || iEnergyBin  == 8 || iEnergyBin  == 9)
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                    }
            else if(iEtaBin == 5)
            {

                        if( iEnergyBin == 0 )
                        {
                            //if(ibias == "9" || ibias == "11" || ibias == "14" || ibias == "15" )
                            //{
                                teste = true;
                                f2->SetParameters( 2.36020e+00  , 1.04245e+00  , 5.98687e-02 ,  -6.89685e-01  , 1.87169e+00);
                            //}else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                        }
                        else if (iEnergyBin == 5 )
                        {
                                //if(iBias == 1 || iBias == 2 || iBias == 8 )
                                if(iBias == 1 || iBias == 2 || iBias == 10 )
                                {
                                    teste = true;
                                    f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                                }
                                else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                        }
                        else if (iEnergyBin == 6 )
                        {
                            if(iBias == 8 || iBias == 10 || iBias == 12)
                                {
                                    teste = true;
                                    f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                                }
                                else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
                        }
                        else if (iEnergyBin == 7 )
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                        }
                        else if (iEnergyBin == 8 )
                        {
                            if(iBias == 11)
                            {
                                teste = true;
                                f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                            }
                            else
                            {
                                teste = true;
                                f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                            }

                        }
                        else if( iEnergyBin == 9    )
                        {
                            if(iBias == 1  || iBias == 3 || iBias == 14 )
                            {
                                teste = true;
                                f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                            }
                            else if(iBias == 10 || iBias == 12 || iBias == 13)
                            {
                                teste = true;
                                f2->SetParameters( 5.28461e-02 ,  1.03349e+00   , 8.38022e-02   , -4.22860e-01 ,  4.38926e+00);
                            }
                            else
                            {
                                teste = true;
                                f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                            }

                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                    }
            else if(iEtaBin == 6)
            {

                        if( iEnergyBin == 0 )
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);

                        }
                        else if (iEnergyBin == 3 && iBias == 5 )
                        {
                                teste = true;
                                f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);

                        }
                        else if (iEnergyBin == 3 && iBias == 7 )
                        {
                                teste = true;
                                f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);

                        }
                        else if (iEnergyBin == 4 )
                        {
                                if(iBias == 8 || iBias == 10 || iBias == 15 || iBias == 17)
                                {
                                    teste = true;
                                    f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                                }else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);


                        }
                        else if (iEnergyBin == 5 && iBias == 1 )
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);

                        }
                        else if (iEnergyBin == 5 && iBias == 3 )
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);

                        }
                        else if (iEnergyBin == 6 )
                        {

                            if( iBias == 0 || iBias == 2 || iBias == 6 || iBias == 8 )
                            {
                                    teste = true;
                                    f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                            }
                            else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                        }
                        else if( iEnergyBin == 7 )
                        {
                            if(iBias == 2 || iBias == 4 || iBias == 6 || iBias == 8 || iBias == 10)
                            {
                                teste = true;
                                f2->SetParameters( 4.91392e-01 ,  1.02770e+00   , 8.45974e-02   , -4.17898e-01  ,  3.99852e+00   );
                            }
                            else
                            {
                                teste = true;
                                f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                            }
                        }
                        else if( iEnergyBin == 8 )
                        {
                            if(iBias == 1 || iBias == 3 || iBias == 8 || iBias == 10 || iBias == 12 || iBias == 13 )
                            {
                                teste = true;
                                f2->SetParameters( 5.28461e-02 ,  1.03349e+00   , 8.38022e-02   , -4.22860e-01 ,  4.38926e+00);
                            }
                            else
                            {
                                teste = true;
                                f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);
                            }
                        }
                        else if( iEnergyBin == 9 )
                        {
                            teste = true;
                            f2->SetParameters( 8.23483e+00  , 9.99187e-01  , 5.03072e-02   ,  -5.25714e-01 , 1.79204e+00);

                        }
                        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);

                    }
            else if(iEtaBin == 7)
            {
                        if(iEnergyBin == 9 && iBias == 4)
                        {
                            teste = true;
                            f2->SetParameters( 1.74703e+01  ,  9.99404e-01 , 1.46469e-01 , -3.24588e-01 , 1.03066e+00 );
                        }
                        else
                        {
                            teste = true;
                            f2->SetParameters( 1.11068e+01  ,  1.04184e+00  , 1.00382e-01   ,  -5.77218e-01 , 4.34526e+00);
                        }

                    }
            */
            else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
        }
        else
        {
             
            if(map_histog->GetEntries() >= 100e3)
            {
                cout<<" >=100e3 " << endl;
                cout<<" entries like " << map_histog->GetEntries() << endl;

                if(isNomFitRange)
                {
                    teste = true;
                    f2->SetParameters( 5.73926e+04 , 9.89075e-01  , 6.19108e-02  ,  -1.03298e+00  ,   1.59785e+00 );
                }
                else
                {
                    teste = true;
                    f2->SetParameters( 2.79598e+03, 9.92913e-01   , 6.39535e-02    , -7.04606e-01,  1.53934e+00 );
                }
                
            }
            else
            {
                cout<<" <=100e3 " << endl;
                cout<<" entries like " << map_histog->GetEntries() << endl;

                if(isNomFitRange)
                {
                        f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00); 
                    }
                else
                {
                        teste = true;
                        f2->SetParameters( 2.79598e+03, 9.92913e-01   , 6.39535e-02    , -7.04606e-01,  1.53934e+00 );
                    }
            
            } 
        }
    
    }
    else if(iperiod == 1)
    {
    cout<<" Entries in period 1 " << map_histog->GetEntries() << endl;
    if(isNomFitRange)
    {
        if(isNomElectronID)
        {
        
            cout<< checkedFailedFit(iEtaBin, iEnergyBin, iBias) << endl;
            bool isFitFailed = checkedFailedFit(iEtaBin, iEnergyBin, iBias);
                
            //bool isFitFailed = false;
            if(isFitFailed)
            {
            
                vector<double> new_params = fixingOnBadFits( iEtaBin, iEnergyBin, iBias );   
                teste = true;
                f2->SetParameters( new_params.at(0) , new_params.at(1)  , new_params.at(2)  ,  new_params.at(3)  ,   new_params.at(4) );
            }
            else
            {
                f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
            }
        }
        else 
        {
            cout<<"electronID"<<endl;
            teste = true;
            f2->SetParameters( 1.00000e+05 ,  1.01322e+00   , 5.28402e-02   ,  -5.82989e-01  , 3.31382e+00  );
        }
        
          
    }
    else
    {
        
        // for 0.95-1.5 only
        if ( iEtaBin == 0 )
        {
            if(iEnergyBin == 3)
            {
                teste = true;
                //f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
                f2->SetParameters( 1.75105e+01  , 1.02008e+00 ,  5.66517e-02  ,  -6.74312e-01   ,  1.52587e+00  );
            }
            else  if(iEnergyBin == 5)
            {
                teste = true;
                f2->SetParameters( 1.75105e+01  , 1.02008e+00 ,  5.66517e-02  ,  -6.74312e-01   ,  1.52587e+00  );
            }
            else if(iEnergyBin == 7 || iEnergyBin == 8 || iEnergyBin == 9)
            {
                             teste = true;
                             f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
            }
            else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
        }
        else if ( iEtaBin == 1 && iEnergyBin == 5 )
        {
            teste = true;
            f2->SetParameters( 9.61259e+04 ,  1.01261e+00   , 4.59479e-02  , -5.92626e-01     ,  1.95821e+00 );
        }
        else if ( iEtaBin == 1 && iEnergyBin == 8 )
        {
            teste = true;
            f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
        }
        else if ( iEtaBin == 2 && iEnergyBin == 6  )
        {
            teste = true;
            f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
        }
        else if ( iEtaBin == 2 && iEnergyBin == 9  )
        {
            teste = true;
            f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
        }
        else if ( iEtaBin == 4 && iEnergyBin == 6  )
        {
                         teste = true;
                         f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
                     }
        else if (iEtaBin == 5)
        {
            if(iEnergyBin == 0)
            {
                teste = true;
                f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 2.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
            }
            else if(iEnergyBin == 8 || iEnergyBin == 9)
            {
                teste = true;
                //f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 5.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
                f2->SetParameters( 5.81249e+01  ,  1.00872e+00  ,2.44957e-02 , -6.75494e-01  , 2.24088e+00);
            }
            else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
        }
        else if (iEtaBin == 6)
        {
          if(iEnergyBin == 0 ) 
           {
                teste = true;
                f2->SetParameters( 5.81249e+01  ,  1.00872e+00  ,6.44957e-02 , -6.75494e-01  , 2.24088e+00);
           }
            else if (iEnergyBin == 3)
            {
                teste = true;
                f2->SetParameters( 9.22727e+04   , 1.01995e+00  , 4.36835e-02   ,  -6.59625e-01   ,  1.81919e+00   );
            }
            else if (iEnergyBin == 7)
            {
                teste = true;
                f2->SetParameters( 5.81249e+01  ,  1.00872e+00  ,6.44957e-02 , -6.75494e-01  , 2.24088e+00);
            }
            else if (iEnergyBin == 8 || iEnergyBin == 9)
            {
                teste = true;
                //f2->SetParameters( 5.81249e+01  ,  1.00872e+00  ,6.44957e-02 , -6.75494e-01  , 2.24088e+00);
                f2->SetParameters( 5.81249e+01  ,  1.00872e+00  ,2.44957e-02 , -6.75494e-01  , 2.24088e+00);
            }
            else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00);
        }
        else if (iEtaBin == 7 && iEnergyBin == 1 )
        {
            teste = true;
            f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 2.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
        }
        else if (iEtaBin == 7 && iEnergyBin == 4 )
        {
            teste = true;
            f2->SetParameters( 3.10914e+00 ,  1.00462e+00   , 2.75953e-02  ,    -5.46311e-01    ,   1.59785e+00 );
        }
        else f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00); //fit range 0.9-1.3
        
        //f2->SetParameters( 2.85935e+02  ,  1.00161e+00    , 4.01616e-02 ,   -8.81297e-01   ,  2.44256e+00); //fit range 0.9-1.3
    }
    
}
    
    if(!noParms)
    {
        //using itif(!teste) f2->SetParLimits(0,1.,1e+5); //nominal fit
        f2->SetParLimits(0,1.,1e+10);
        f2->SetParLimits(1,0.8,1.2);
        f2->SetParLimits(2,0.01,0.5);
        f2->SetParLimits(3,-2.,2.);
        f2->SetParLimits(4,1.,10);
    }
    
    TH1F *hRatio = new TH1F(Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i",iperiod,iEtaBin,iEnergyBin,iBias), "", 50, 0.6, 2.); 
    hRatio->Sumw2();

    gStyle->SetOptFit(1);
    gROOT->ForceStyle();

    TCanvas* canvas = new TCanvas("canvas","canvas",800,600) ;
    canvas->cd() ; 
    TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
    TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
    pad1->Draw();pad2->Draw();
    pad1->cd();
    pad1->SetBottomMargin(0.16);
    pad2->SetBottomMargin(0.24);
    map_histog->SetMarkerStyle(20);
    map_histog->SetMarkerSize(1.);
    map_histog->GetYaxis()->SetTitle("Events");
    map_histog->GetXaxis()->SetTitle("EoverP");
    //map_histog->UseCurrentStyle();
    map_histog->Draw("E");
    map_histog->SaveAs(path+Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i.root",iperiod,iEtaBin,iEnergyBin,iBias));
    map_histog->Fit("f2", "r");

    //map_histog->Print("all");
    
    

    //if(iEtaBin <= 7){
    if(iEtaBin <= 16){
    cout<< " low ETA RATIO " << endl;
    //map_histog->Print("all");
        
    double sumDelta = 0;
    double newRatio = 0;
        
    for (int ibin = 11; ibin <= 25; ++ibin) {// nominal fit range
    int low, high;
    //low = 13; high = 32; 0.95-1.5
    /*if(iEtaBin < 4)
    {
        low = 13; high = 32;
    }else
    {
        low = 10; high = 32;
    }*/
    //for (int ibin = low; ibin <= high; ++ibin) {
        if( map_histog->GetBinContent(ibin) > 0 ){
                    double res =  (map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )/map_histog->GetBinError(ibin);
                    hRatio->SetBinContent(ibin, res  );
                    hRatio->SetBinError(ibin, map_histog->GetBinError(ibin)  );
                    //cout<<" data value  " << map_histog->GetBinContent(ibin) << " function value " << f2->Eval( map_histog->GetBinCenter(ibin) ) << " diff " <<  fabs((map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )) << endl;
                    sumDelta += (fabs((map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )));
            
            
        }
   }   
        
   newRatio = sumDelta/((25-11)+1);
   cout<<" NEW RATIO " << newRatio << endl;
        
   TString status_fit = gMinuit->fCstatu;

   file_brokenfit << iEtaBin <<"\t" << iEnergyBin << "\t" <<iBias << "\t" << newRatio << "\t"<<f2->GetParameter(0) << "\t" << f2->GetParameter(1) << "\t" << f2->GetParameter(2) << "\t" << f2->GetParameter(3) << "\t" << f2->GetParameter(4) << endl;
   
        
        
} 
    else{
    cout<< " higher ETA RATIO " << endl;
    for (int ibin = 10; ibin <= 32; ++ibin) {
        if( map_histog->GetBinContent(ibin) > 0 ){
                    double res =  (map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )/map_histog->GetBinError(ibin);
                    hRatio->SetBinContent(ibin, res  );
                    hRatio->SetBinError(ibin, map_histog->GetBinError(ibin)  );
        }
   }       
}
    
    pad2->cd();
  hRatio->GetXaxis()->SetTitle((Form("EoverP, #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));
    
  //hRatio->GetXaxis()->SetTitle((Form("EoverP, #eta bin %i: abs(#eta) in range [%g, %g], E1/E2", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1] )));    
    
  hRatio->GetYaxis()->SetTitle("DATA / FIT");
  //hRatio->GetXaxis()->SetRangeUser(84., 98.);
  //hRatio->GetYaxis()->SetRangeUser(0.,1.3);
  hRatio->GetYaxis()->SetRangeUser(-50.,50);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.08);
  hRatio->GetXaxis()->SetTitleSize(0.08);
  hRatio->GetYaxis()->SetTitleSize(0.09);
  hRatio->GetYaxis()->SetTitleOffset(0.6);
  hRatio->GetXaxis()->SetTitleOffset(1.2);
  hRatio->Draw("E1");
  //pad2->Modified();
  //pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();   
  
  canvas->SaveAs(path+Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i.pdf",iperiod,iEtaBin,iEnergyBin,iBias));

 
      delete pad1;
      delete pad2;
      delete canvas;

     double meanValue, errMeanValue;

    if(iperiod == 0)
    {
        if(iEtaBin == 12 && iEnergyBin > 6 )
        {
            meanValue = 0;
            errMeanValue  = 0;
        }
        else
        {
            meanValue = f2->GetParameter(1);
            errMeanValue = f2->GetParError(1);
        }
    }
    else if(iperiod == 1)
    {
        meanValue = f2->GetParameter(1);
        errMeanValue = f2->GetParError(1);
    }

     //double meanValue = f2->GetParameter(1);
     //double errMeanValue = f2->GetParError(1);

     parameters.push_back(meanValue);//0
     parameters.push_back(errMeanValue);//1


    //cout<<" meanValue " << meanValue  << " errMeanValue  " << errMeanValue  << endl;

    //return make_pair(meanValue, errMeanValue);
}
else
{
    static bool bInitialized = false;
    if(!bInitialized)
    {
        Initialize();
        bInitialized = true;
    }
    
    vector<float> param;
    vector<std::shared_ptr<FitInfo>> v;
    if(iBias != -1)
    {
        v = GetFitInfo(iEtaBin, iEnergyBin, true);
    }
    else
    {
        if(iperiod == 0)
            v = GetFitInfoExperimentalData(iEtaBin);
        if(iperiod == 1)
            v = GetFitInfoMCSimulation(iEtaBin);
    }
    
    for(int i = 0; i < v.size(); i++)
    {
        if(v[i]->broken && v[i]->bias == iBias)
        {
            cout << "bin: " << v[i]->eta << " energy: " << v[i]->energy << " bias: " << v[i]->bias << " period: " << iperiod << " is broken" << endl;
            //it needs to be fixed
            std::vector<std::shared_ptr<LinearFitFunc>> funcs;
            if(iBias != -1 && gFixFuncMap.find(std::make_pair(iEtaBin, iEnergyBin)) != gFixFuncMap.end())
                funcs = gFixFuncMap[std::make_pair(iEtaBin, iEnergyBin)];
            else if(iperiod == 0 && gFixFuncDataMap.find(iEtaBin) != gFixFuncDataMap.end())
                funcs = gFixFuncDataMap[iEtaBin];
            else if(iperiod == 1 && gFixFuncMCMap.find(iEtaBin) != gFixFuncMCMap.end())
                funcs = gFixFuncMCMap[iEtaBin];
            else
                ;
            
            if(!funcs.empty())
            {
                double x = (iBias == -1) ? iEnergyBin : iBias;
                //we found a fix function
                param.push_back(funcs[0]->Func(x));
                param.push_back(funcs[1]->Func(x));
                param.push_back(funcs[2]->Func(x));
                param.push_back(funcs[3]->Func(x));
                param.push_back(funcs[4]->Func(x));
                cout << "fixed value is :" << param[0] << " " << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << std::endl;
                           
                f2->SetParameters( param[0] , param[1]  , param[2] ,  param[3]  , param[4] );
                           
                f2->SetParLimits(0,1.,1e+10);
                f2->SetParLimits(1,0.8,1.2);
                f2->SetParLimits(2,0.01,0.5);
                f2->SetParLimits(3,-2.,2.);
                f2->SetParLimits(4,1.,10);
                           
                TH1F *hRatio = new TH1F(Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i",iperiod,iEtaBin,iEnergyBin,iBias), "", 50, 0.6, 2.);
                hRatio->Sumw2();

                gStyle->SetOptFit(1);
                gROOT->ForceStyle();

                TCanvas* canvas = new TCanvas("canvas","canvas",800,600) ;
                canvas->cd() ;
                TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
                TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
                pad1->Draw();pad2->Draw();
                pad1->cd();
                pad1->SetBottomMargin(0.16);
                pad2->SetBottomMargin(0.24);
                map_histog->SetMarkerStyle(20);
                map_histog->SetMarkerSize(1.);
                map_histog->GetYaxis()->SetTitle("Events");
                map_histog->GetXaxis()->SetTitle("EoverP");
                //map_histog->UseCurrentStyle();
                map_histog->Draw("E");
                map_histog->SaveAs(path+Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i.root",iperiod,iEtaBin,iEnergyBin,iBias));
                map_histog->Fit("f2", "r");

                TString status_fit = gMinuit->fCstatu;

                file_brokenfit << "eta" << iEtaBin <<"\t"<<"energy" << iEnergyBin <<"\t"<< "bias"<<iBias << "\t"<< status_fit << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << " " << f2->GetParameter(3) << "  " << f2->GetParameter(4) << endl;

                //if(iEtaBin <= 7){
                if(iEtaBin <= 16)
                {
                    cout<< " low ETA RATIO " << endl;
                    //map_histog->Print("all");
                    for (int ibin = 11; ibin <= 32; ++ibin)
                    {
                        // nominal fit range
                        int low, high;
                        //low = 13; high = 32; 0.95-1.5
                        //for (int ibin = low; ibin <= high; ++ibin) {
                        if( map_histog->GetBinContent(ibin) > 0 )
                        {
                            double res =  (map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )/map_histog->GetBinError(ibin);
                            hRatio->SetBinContent(ibin, res  );
                            hRatio->SetBinError(ibin, map_histog->GetBinError(ibin)  );
                        }
                    }
                }
                else
                {
                    cout<< " higher ETA RATIO " << endl;
                    for (int ibin = 10; ibin <= 32; ++ibin)
                    {
                        if( map_histog->GetBinContent(ibin) > 0 )
                        {
                            double res =  (map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )/map_histog->GetBinError(ibin);
                            hRatio->SetBinContent(ibin, res  );
                            hRatio->SetBinError(ibin, map_histog->GetBinError(ibin)  );
                        }
                    }
                    
                }
                   
                pad2->cd();
                hRatio->GetXaxis()->SetTitle((Form("EoverP, #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));

                hRatio->GetYaxis()->SetTitle("DATA / FIT");
                hRatio->GetYaxis()->SetRangeUser(-50.,50);
                hRatio->GetXaxis()->SetLabelSize(0.1);
                hRatio->GetYaxis()->SetLabelSize(0.08);
                hRatio->GetXaxis()->SetTitleSize(0.08);
                hRatio->GetYaxis()->SetTitleSize(0.09);
                hRatio->GetYaxis()->SetTitleOffset(0.6);
                hRatio->GetXaxis()->SetTitleOffset(1.2);
                hRatio->Draw("E1");
                pad2->RedrawAxis();
                pad2->Update();

                canvas->SaveAs(path+Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i.pdf",iperiod,iEtaBin,iEnergyBin,iBias));

                
                delete pad1;
                delete pad2;
                delete canvas;

                double meanValue = f2->GetParameter(1);
                double errMeanValue = f2->GetParError(1);

                parameters.push_back(meanValue);//0
                parameters.push_back(errMeanValue);//1
                return parameters;
            }
            else
            {
                cout << "we don't find an availabe fix function for eta: " << iEtaBin << " energy: " << iEnergyBin << " bias: " << iBias << " period: " << iperiod << endl;
            }
        }
    }
    
    parameters.push_back(0);//0
    parameters.push_back(0);//1    
    return parameters;
}
     
return parameters;

 }
//////////////////////////////////////////SAVING HISTOS////////////////////////////////////////





void mySave(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, int iEnergyBin,TString vartagname, TString varname,TString path, int m_layer){	

  gStyle->SetOptFit(0);
  TCanvas *canvas = new TCanvas("scanning","scanning", 800,800);
  canvas->Clear();
  canvas->cd();
  TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
  TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
  pad1->Draw();pad2->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.16);
  pad2->SetBottomMargin(0.24);
  //gPad->SetLogy();

  
  map_histog_data->GetYaxis()->SetTitle("Normalised Events");
  map_histog_data->GetXaxis()->SetRangeUser(lowcut, highcut);
  map_histog_data->GetYaxis()->SetRangeUser(-0.1, 2*map_histog_data->GetMaximum());
  map_histog_data->GetYaxis()->SetTitleSize(0.07);
  map_histog_data->GetYaxis()->SetTitleOffset(0.9);
  map_histog_data->DrawNormalized("E1");
  map_histog_mc->SetLineColor(2);
  map_histog_mc->SetMarkerColor(2);
  map_histog_mc->DrawNormalized("histsame");
  
  TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader(Form("#eta bin %i: abs(#eta) in range [%g, %g]", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1]));
  

  leg->AddEntry(map_histog_mc,"MC","pl");
  if(iBias < 2)
    leg->AddEntry(map_histog_data, "DATA","l");
  else
    leg->AddEntry(map_histog_data, Form("DATA_bias%i",iBias-2),"l");

  pad1->Modified();
  pad1->RedrawAxis();
  pad1->Update();
  pad2->cd();
  
  TH1F* clone_mc = (TH1F*)map_histog_mc;
  TH1F* clone_data = (TH1F*)map_histog_data;
  clone_mc->Scale(clone_data->Integral()/clone_mc->Integral());
  clone_data->Divide(clone_mc);
  
  clone_data->GetXaxis()->SetTitle(vartagname);
  clone_data->GetYaxis()->SetTitle("DATA / MC");
  clone_data->GetXaxis()->SetRangeUser(lowcut, highcut);
  clone_data->GetYaxis()->SetRangeUser(0.7,1.3);
  clone_data->GetXaxis()->SetLabelSize(0.1);
  clone_data->GetYaxis()->SetLabelSize(0.08);
  clone_data->GetXaxis()->SetTitleSize(0.08);
  clone_data->GetYaxis()->SetTitleSize(0.09);
  clone_data->GetYaxis()->SetTitleOffset(0.6);
  clone_data->GetXaxis()->SetTitleOffset(1.2);
  clone_data->Draw("E1");
  pad2->Modified();
  pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();

  TString output_folder = path+Form("calibrated_%s",varname.Data());
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } 
  else {
    gSystem->mkdir(path+Form("calibrated_%s/",varname.Data()), kTRUE);
  }

  if (!(stat(output_folder+Form("/Histograms/layer%i/",m_layer), &sb) == 0 && S_ISDIR(sb.st_mode))) 
    gSystem->mkdir(output_folder+Form("/Histograms/layer%i/",m_layer), kTRUE);

  if(iBias<2){
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s.pdf",m_layer,iEtaBin, iEnergyBin, varname.Data()));
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s.root",m_layer,iEtaBin, iEnergyBin, varname.Data()));
  }
  else {
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s_bias%i.pdf",m_layer,iEtaBin, iEnergyBin, varname.Data(),iBias-2));
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s_bias%i.root",m_layer,iEtaBin, iEnergyBin, varname.Data(),iBias-2));
  }

  delete canvas;      
  
}



void mySave2D(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, TString vartagname, TString varname, TString path,TString parameter, int m_layer){	


  cout<<" MY SAVE 2D " << endl;    
  /*    
  for(int i = 0; i < map_histog_mc->GetNbinsX(); i++)
  {
      cout<<" bin " << i << " get bin content MC " << map_histog_mc->GetBinContent(i) << " get bin content MC var " << map_histog_data->GetBinContent(i) << endl;
  }
  */
    
  gStyle->SetOptFit(0);
  TCanvas *canvas2D = new TCanvas("scanning","scanning", 800,800);
  canvas2D->Clear();
  canvas2D->cd();
  if(varname=="m_mass2el") map_histog_mc->GetYaxis()->SetTitle("m_{ee} peak [GeV]");
  if(varname=="m_eOverp") map_histog_mc->GetYaxis()->SetTitle("E/p peak");
  map_histog_mc->GetXaxis()->SetTitle("E1/E2");
  map_histog_mc->SetLineColor(2);
  map_histog_mc->SetMarkerColor(2);
  map_histog_data->SetMarkerColor(1);
  map_histog_data->SetLineColor(1);
  if(varname=="m_mass2el") map_histog_mc->GetYaxis()->SetRangeUser(lowcut, highcut);
  if(varname=="m_eOverp")  map_histog_mc->GetYaxis()->SetRangeUser(0., 5.);

  map_histog_mc->Draw("E1");
  map_histog_data->Draw("E1 samese");
  TLegend *leg2D = new TLegend(0.5,0.6,0.9,0.9);
  if(varname=="m_mass2el") leg2D->SetHeader(Form("#eta bin %i: abs(#eta) in range [%g, %g]", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1]));
  leg2D->SetLineColor(0);
  leg2D->SetShadowColor(0);
  leg2D->SetTextSize(0.035);
  leg2D->AddEntry(map_histog_data,"DATA");
  leg2D->AddEntry(map_histog_mc, "MC");
  leg2D->Draw("");


  TString output_folder = path+Form("calibrated_%s",varname.Data());
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } else {
    gSystem->mkdir(path+Form("calibrated_%s/",varname.Data()), kTRUE);
  }

  if (!(stat(output_folder+Form("/Histograms2D/layer%i/",m_layer), &sb) == 0 && S_ISDIR(sb.st_mode))) 
    gSystem->mkdir(output_folder+Form("/Histograms2D/layer%i/",m_layer), kTRUE);


  if(iBias<2){  
    map_histog_mc->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo%s_mc_%i_%s.root",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    map_histog_data->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo%s_data_%i_%s.root",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_%s.pdf",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_%s.root",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    
    //save the parameters for MC
    if(varname=="m_mass2el" || varname=="m_eOverp"){      
      map_histog_mc->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_mc"+parameter+"_%i_%s.root",m_layer,iEtaBin,varname.Data()));
      map_histog_data->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_data"+parameter+"_%i_%s.root",m_layer,iEtaBin,varname.Data()));
    }
  }
	    
  else{  
      
    cout<<"FABIOLIFE"<<endl; 
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_bias%i_%s.pdf",m_layer,parameter.Data(),iEtaBin,iBias,varname.Data()));
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_bias%i_%s.root",m_layer,parameter.Data(),iEtaBin,iBias,varname.Data()));
    
    //save the parameters for MC
    if(varname=="m_mass2el" || varname=="m_eOverp"){      
      map_histog_mc->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_mc"+parameter+"_bias%i_%i_%s.root",m_layer,iEtaBin,iBias,varname.Data()));
      map_histog_data->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_data"+parameter+"_bias%i_%i_%s.root",m_layer,iEtaBin,iBias,varname.Data()));
    }
  }  
}



vector<double> fixingOnBadFits( int iEtaBin, int iEnergyBin, int iBias )
{

    
    TString path_fitRes = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/";
    
    ifstream file( path_fitRes+"out_fit2.dat"  );
    
    vector<TString> coordinates;
    vector<TString> fit_conv;
    vector<double> parameters;
    
    
    
    
    std::string str; 
    std::string foundConverge;
    
    bool found = false;
    
    while (std::getline(file, str))
    {
        TString findFailed = str;
        TString foundConv;

        if( findFailed.Contains("FAILED") || findFailed.Contains("NOT POSDEF") )
        {
                TString eta = Form("eta%d",iEtaBin);
                TString energy = Form("energy%d",iEnergyBin);
                TString bias = Form("bias%d",iBias);
                //cout<<" ID failed " << eta << "\t" << energy << "\t" << bias << endl;
            
            

                std::string tmp_str;    
                while (std::getline(file, tmp_str))
                {
                    
                    TString findConv = tmp_str; 
                    if( findConv.Contains(eta) && findConv.Contains(energy) )
                    {
                        
                        if(findConv.Contains("CONVERGED"))
                        {   
                            cout<<" CONVERGE FOUND " << findConv << endl; 
                            foundConverge = tmp_str; 
                            break;
                            
                        }
                    }
                    break;
                
                }
                
            }
        
            
    }
    
  
    std::istringstream ss(foundConverge);
    std::vector<string> vectors;
    std::vector<double> params;
                    
    string v;
    while (ss)                
    { 
        ss >> v;
        vectors.push_back(v);
    }
                          
    
    cout<<" Let's use: " << foundConverge << "\t" << vectors.size() << endl;
    
    int i = 4;
    while(i < 9)
    {
        double convertedValue = std::stod(vectors.at(i));
        params.push_back(convertedValue);
        i++;
    }
    
    cout<<" SIZE " << params.size() << endl;
    return params;
}

bool checkedFailedFit( int iEtaBin, int iEnergyBin, int iBias )
{

    
    TString path_fitRes = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/";
    ifstream file( path_fitRes+"out_fit_2.dat"  );
    
    
    std::string str; 
    std::string foundConverge;
    
    bool fitStatus = false;
    
    while (std::getline(file, str))
    {
        TString findFailed = str;
        if(findFailed.Contains("FAILED") || findFailed.Contains("NOT POSDEF") )
        {
            fitStatus = true;
            break;
                
        }   
    }
    cout<<" fitStatus " << fitStatus << endl;
    
  
    
    return fitStatus;
}

void simpleOnBadFits( int iEtaBin, int iEnergyBin, int iBias )
{

    
    TString path_fitRes = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/";
    
    ifstream file( path_fitRes+"out_fit2.dat"  );
    ifstream file2( path_fitRes+"out_fit2.dat"  );
    
    
    TString eta = Form("eta%d",iEtaBin);
    TString energy = Form("energy%d",iEnergyBin);
    TString bias = Form("bias%d",iBias);
    
    TString newbias = Form("bias%d",iBias-1);

    std::string str; 
    std::string str2; 
    
    double constParam;
    double constParam2;
    //double constParam;
    while (std::getline(file, str))
    {
        TString lineTobeFound = str;
        
        if( lineTobeFound.Contains(eta) && lineTobeFound.Contains(energy) && lineTobeFound.Contains(bias) )
        {
            cout << " lineTobeFound " << lineTobeFound << endl;
           
            TString findFailed = str;
            TString foundConv;

            std::istringstream ss1(str);
            std::vector<string> tmp1;
                        
            string st1;
            while (ss1)                
            { 
                ss1 >> st1;
                tmp1.push_back(st1);
            }
        
            constParam = std::stod(tmp1.at(4));   
        }
            
            
        
    
    }
    
    while (std::getline(file, str))
    {
        cout<<" FABIO " << endl;
        std::istringstream ss1(str);
        std::vector<string> tmp1;
                        
            string st1;
            while (ss1)                
            { 
                ss1 >> st1;
                tmp1.push_back(st1);
            }
        
            constParam2 = std::stod(tmp1.at(4));  
            cout<<" constParam2 " << constParam2 << endl;
    }
    
   
}

void LoadFitInfo(string path)
{
    if(gInfo.empty())
    {
        cout<<" Read fitinto from:  " << path<< endl;
        ifstream file2(path);
        std::string str;
        int bin, energy, bias;
        double error, constant, mean, sigma, alpha, n;
        vector<FitInfo> res;
        
        if(!file2.is_open())
        {
            cout << "the state of the filestream is not right" <<endl;
            cout << "Error: " << strerror(errno) <<endl;
            return;
        }
        
        
        while (std::getline(file2, str))
        {
            std::istringstream ss(str);
            std::vector<std::string> numbers((std::istream_iterator<std::string>(ss)),
                                     std::istream_iterator<std::string>());
            bin = stoi(numbers[0]);
            energy = stoi(numbers[1]);
            bias = stoi(numbers[2]);
            error = stod(numbers[3]);
            constant = stod(numbers[4]);
            mean = stod(numbers[5]);
            sigma = stod(numbers[6]);
            alpha = stod(numbers[7]);
            n = stod(numbers[8]);
                
            gInfo.push_back(std::make_shared<FitInfo>(bin, energy, bias, constant, mean, sigma, alpha, n, error));
        }
        cout<< "Read in total: " << gInfo.size() << "fit info data" << endl;
        file2.close();
    }
}

//Get all the fitinfo belong to iEtaBin and iEnergyBin, but we don't include FitInfo with bias = -1
vector<std::shared_ptr<FitInfo>> GetFitInfo(int iEtaBin, int iEnergyzBin, bool excludingBiasMinusOne)
{
    vector<std::shared_ptr<FitInfo>> res;
    //load the fitinfo if it is not loaded yet
    LoadFitInfo("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/dataPoints.dat");
    for(int i = 0; i < gInfo.size(); i++)
    {
        if(gInfo[i]->eta == iEtaBin && gInfo[i]->energy == iEnergyzBin)
        {
            if(gInfo[i]->bias == -1 && excludingBiasMinusOne)
                continue;
            res.push_back(gInfo[i]);
        }
    }
    return res;
}

//only get the simulation data with eta = iEtaBin and bias = -1
vector<std::shared_ptr<FitInfo>> GetFitInfoMCSimulation(int iEtaBin)
{
    vector<std::shared_ptr<FitInfo>> res;
    LoadFitInfo("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/dataPoints.dat");
    double read = true;
    for(int i = 0; i < gInfo.size(); i++)
    {
        if(gInfo[i]->eta == iEtaBin)
        {
            if(gInfo[i]->bias == -1)
            {
                if(read)
                {
                    res.push_back(gInfo[i]);
                    read = false;
                }
                else
                {
                    read = true;
                }
            }
        }
    }
    return res;
}

vector<std::shared_ptr<FitInfo>> GetFitInfoExperimentalData(int iEtaBin)
{
    vector<std::shared_ptr<FitInfo>> res;
    LoadFitInfo("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/dataPoints.dat");
    double read = false;
    for(int i = 0; i < gInfo.size(); i++)
    {
        if(gInfo[i]->eta == iEtaBin)
        {
            if(gInfo[i]->bias == -1)
            {
                if(read)
                {
                    res.push_back(gInfo[i]);
                    read = false;
                }
                else
                {
                    read = true;
                }
            }
        }
    }
    return res;
}

//this function load the fitinfo, mark broken data and generate fix function
void Initialize()
{
    cout << "begin to initialize" << endl;
    LoadFitInfo("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/dataPoints.dat");
    
    for(int bin = 0; bin <= 15; bin++)
    {
        for(int energy = 0; energy <= 9; energy++)
        {
            vector<std::shared_ptr<FitInfo>> v = GetFitInfo(bin, energy, true);
            MarkBrokenFit(v,  2.5);
            vector<std::shared_ptr<LinearFitFunc>> funcs = GenerateFitFunc(v);
            if(!funcs.empty())
            {
                gFixFuncMap[std::make_pair(bin, energy)] = funcs;
            }
        }
    }
    
    for(int bin = 0; bin <= 15; bin++)
    {
        vector<std::shared_ptr<FitInfo>> v = GetFitInfoExperimentalData(bin);
        MarkBrokenFit(v,  2.5);
        vector<std::shared_ptr<LinearFitFunc>> funcs = GenerateFitFuncBiasMinusOne(v);
        if(!funcs.empty())
        {
            gFixFuncDataMap.insert(std::make_pair(bin, funcs));
        }
    }
    
    for(int bin = 0; bin <= 15; bin++)
    {
        vector<std::shared_ptr<FitInfo>> v = GetFitInfoMCSimulation(bin);
        MarkBrokenFit(v,  2.5);
        vector<std::shared_ptr<LinearFitFunc>> funcs = GenerateFitFuncBiasMinusOne(v);
        if(!funcs.empty())
        {
            gFixFuncMCMap.insert(std::make_pair(bin, funcs));
        }
    }
    cout << "end of initialization" << endl;
}

vector<std::shared_ptr<LinearFitFunc>> GenerateFitFuncBiasMinusOne(const vector<std::shared_ptr<FitInfo>>& info)
{
    vector<std::shared_ptr<LinearFitFunc>> res;
    int numGood = 0;
    for(int i = 0; i < info.size(); i++)
    {
        if(!info[i]->broken)
            numGood++;
    }
    
    if(numGood < 2)
    {
        cout << "Failed to find enough good data points to generate Fit function" << endl;
        return res;
    }
    
    if(numGood == info.size())
    {
        cout << "All are good, no need to generate fit function" << endl;
        return res;
    }
    
    double accumEnergy = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumEnergy += (a->broken ? 0.0 : a->energy);
    });
    double energyMean = accumEnergy / numGood;
    
    // compute constant fit function
    double accumConstant = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumConstant+= (a->broken ? 0.0 : a->constant);
    });
    double constantMean = accumConstant / numGood;
    
    double accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->energy - energyMean) * (a->constant - constantMean);
    });
    
    double accum2 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum2 += a->broken ? 0.0 : (a->energy - energyMean) * (a->energy - energyMean);
    });
    
    double betaConstant = accum1 / accum2;
    double alphaConstant = constantMean - betaConstant * energyMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaConstant, betaConstant));
    //compute mean fit func
    double accumMean = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumMean+= (a->broken ? 0.0 : a->mean);
    });
    double meanMean = accumMean / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->energy - energyMean) * (a->mean - meanMean);
    });
    
    double betaMean = accum1 / accum2;
    double alphaMean = meanMean - betaMean * energyMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaMean, betaMean));
    //compute sigma fit func
    double accumSigma = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumSigma+= (a->broken ? 0.0 : a->sigma);
    });
    double sigmaMean = accumSigma / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->energy - energyMean) * (a->sigma - sigmaMean);
    });
    
    double betaSigma = accum1 / accum2;
    double alphaSigma = sigmaMean - betaSigma * energyMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaSigma, betaSigma));
    
    //compute alpha fit func
    double accumAlpha = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumAlpha += (a->broken ? 0.0 : a->alpha);
    });
    double alphaMean1 = accumAlpha / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->energy - energyMean) * (a->alpha - alphaMean1);
    });
    
    double betaAlpha = accum1 / accum2;
    double alphaAlpha = alphaMean1 - betaAlpha * energyMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaAlpha, betaAlpha));
    
    //compute n fit func
    double accumN = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumN += (a->broken ? 0.0 : a->n);
    });
    double nMean = accumN / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->energy - energyMean) * (a->n - nMean);
    });
    
    double betaN = accum1 / accum2;
    double alphaN = nMean - betaN * energyMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaN, betaN));
    return res;
}

vector<std::shared_ptr<LinearFitFunc>> GenerateFitFunc(const vector<std::shared_ptr<FitInfo>>& info)
{
    vector<std::shared_ptr<LinearFitFunc>> res;
    bool numGood = 0;
    for(int i = 0; i < info.size(); i++)
    {
        if(!info[i]->broken)
            numGood++;
    }
    
    if(numGood < 2)
    {
        cout << "Failed to find enough good data points to generate Fit function" << endl;
        return res;
    }
    
    if(numGood == info.size())
    {
        cout << "All are good, no need to generate fit function" << endl;
        return res;
    }
    
    double accumBias = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumBias += (a->broken ? 0.0 : a->bias);
    });
    double biasMean = accumBias / numGood;
    
    // compute constant fit function
    double accumConstant = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumConstant+= (a->broken ? 0.0 : a->constant);
    });
    double constantMean = accumConstant / numGood;
    
    double accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->bias - biasMean) * (a->constant - constantMean);
    });
    
    double accum2 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum2 += a->broken ? 0.0 : (a->bias - biasMean) * (a->bias - biasMean);
    });
    
    double betaConstant = accum1 / accum2;
    double alphaConstant = constantMean - betaConstant * biasMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaConstant, betaConstant));
    //compute mean fit func
    double accumMean = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumMean+= (a->broken ? 0.0 : a->mean);
    });
    double meanMean = accumMean / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->bias - biasMean) * (a->mean - meanMean);
    });
    
    double betaMean = accum1 / accum2;
    double alphaMean = meanMean - betaMean * biasMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaMean, betaMean));
    //compute sigma fit func
    double accumSigma = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumSigma+= (a->broken ? 0.0 : a->sigma);
    });
    double sigmaMean = accumSigma / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->bias - biasMean) * (a->sigma - sigmaMean);
    });
    
    double betaSigma = accum1 / accum2;
    double alphaSigma = sigmaMean - betaSigma * biasMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaSigma, betaSigma));
    
    //compute alpha fit func
    double accumAlpha = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumAlpha += (a->broken ? 0.0 : a->alpha);
    });
    double alphaMean1 = accumAlpha / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->bias - biasMean) * (a->alpha - alphaMean1);
    });
    
    double betaAlpha = accum1 / accum2;
    double alphaAlpha = alphaMean1 - betaAlpha * biasMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaAlpha, betaAlpha));
    
    //compute n fit func
    double accumN = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accumN += (a->broken ? 0.0 : a->n);
    });
    double nMean = accumN / numGood;
    
    accum1 = 0.0;
    std::for_each (std::begin(info), std::end(info), [&](const std::shared_ptr<FitInfo> a) {
        accum1 += a->broken ? 0.0 : (a->bias - biasMean) * (a->n - nMean);
    });
    
    double betaN = accum1 / accum2;
    double alphaN = nMean - betaN * biasMean;
    res.push_back(std::make_shared<LinearFitFunc>(alphaN, betaN));
    return res;
}

//This is an improved version to detect broken fit using Z-score algorithm
void MarkBrokenFit(vector<std::shared_ptr<FitInfo>>& v,  double threshold)
{
    if(v.size() <= 1)
    {
        cout << "There is no enough data for us to analyze!" << endl;
        return;
    }
    
    //first we compute the mean value of the error
    auto lambda = [&](const std::shared_ptr<FitInfo> a, const std::shared_ptr<FitInfo> b){ return a->error + b->error; };
    double sum = std::accumulate(v.begin(), v.end(), 0.0, lambda);
    double mean = sum / v.size();
    
    //now compute the standard deviation
    double accum = 0.0;
    std::for_each (std::begin(v), std::end(v), [&](const std::shared_ptr<FitInfo> info) {
        accum += (info->error - mean) * (info->error - mean);
    });

    double stdev = sqrt(accum / (v.size() - 1));
    //now mark the broken ones
    for(int i = 0; i < v.size(); i++)
    {
        double z_score = fabs((v[i]->error - mean) / stdev);
        if(z_score > threshold)
            v[i]->broken = true;
    }
}


void FindBrokenBias(vector<FitInfo>& v)
{   
    cout<<" try to find broken for bin " << v[0].eta << "energy: " << v[0].energy << endl;
    double range1 = 0.0;
    double range2 = 0.0;
    for(int i = 0; i < v.size(); i++)
    {
        if(v[i].bias == 2)
            range1 = v[i].error;
        else if(i == v.size() - 1)
        {
            range2 = v[i].error;
        }
        else
            ;
    }
    
    double lower_range;
    double bigger_range;
    if(range1 < range2)
    {
        lower_range = range1;
        bigger_range = range2;
    }
    else
    {
        lower_range = range2;
        bigger_range = range1;
    }
    
    cout<<" lowerrange " << lower_range << " higher range: " << bigger_range << endl;
    lower_range -= lower_range * 0.2;
    bigger_range += bigger_range * 0.2;
    
    for(int i = 0; i < v.size(); i++)
    {
        //cout<<" i: " << i <<  " error: " << v[i].error << endl;
        if(!(v[i].error >= lower_range && v[i].error <= bigger_range))
        {
            //cout<<" index: " << i << "is broken" << endl;
            v[i].broken = true;
        }
    }
    
    cout<<" finish check FindBrokenBias " << endl;
    
}

bool throwKinematics(int iEtaBin, int iEnergyBin, int iBias)
    {
        TString eta = Form("eta%d",iEtaBin);
        TString energy = Form("energy%d",iEnergyBin);
        TString bias = Form("bias%d",iBias);
        
        cout<<" etaf " << eta << " energyf " << energy << " biasf " << bias << endl;
    
        ifstream file1( "./fixThoseBins.dat"  ); 
        std::string str; 

        cout<<" teste 1 " << endl;
    
        TString foundLine;
        bool doFix = false;
        if(file1)
        {
                
            cout<<" teste 1111 " << endl;
            
            while (std::getline(file1, str))
        {
            cout<<" teste 2 " << str << endl;
            TString foundLine = str;
            if(foundLine.Contains(eta) && foundLine.Contains(energy) && foundLine.Contains(bias))
            {
                
                doFix = true;
                break;
            }
                
        }
        }
       
          
        return doFix;
    }


/*
vector<FitInfo> ReadFitData(string fileName)
{
    vector<FitInfo> res;
    ifstream fs(fileName);
}
*/
