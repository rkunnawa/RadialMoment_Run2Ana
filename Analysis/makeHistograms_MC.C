//! Raghav Kunnawalkam Elayavalli
//! Jan 12th 2017
//! Rutgers 
//! for questions or comments: raghav.k.e at CERN dot CH

//! Macro to read MC hiforest and make simple histograms 

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"


double ptbins_jec[] = {15., 18.,  21.,  24.,  28.,  32.,  37.,  43.,   49.,   56.,   64.,   74.,   84.,   97.,  114.,  133.,  153.,  174.,  196.,  220.,  245.,  272.,  300.,  330.,  362.,  395.,  430.,  468.,  507.,  548.,  592.,  638.,  686., 1000.};
const int nbins_pt_jec = sizeof(ptbins_jec)/sizeof(double)-1;

const double etabins_jec[] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, 
            -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
            -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
            -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
            -0.435, -0.348, -0.261, -0.174, -0.087, 
            +0.000, 
            +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
            +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
            +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, 
            +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, 
            +4.363, +4.538, +4.716, +4.889, +5.191
};
const int nbins_eta_jec = sizeof(etabins_jec)/sizeof(double) -1;

static const int trigValue = 4;
static const char trigName [trigValue][256] = {"HLT40","HLT60","HLT80","Combined"};

// divide by bin width
void divideBinWidth(TH1 *h){
  h->Sumw2();
  for (int i=0;i<=h->GetNbinsX();i++){
    Float_t val = h->GetBinContent(i);
    Float_t valErr = h->GetBinError(i);
    val/=h->GetBinWidth(i);
    valErr/=h->GetBinWidth(i);
    h->SetBinContent(i,val);
    h->SetBinError(i,valErr);
  }//binsX loop 
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

float deltaphi(float phi1, float phi2)
{
  float pi=TMath::Pi();
 
  float dphi = TMath::Abs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;

  return TMath::Abs(dphi);
}

double pthatWeight(int startfile = 0){
  // no of events in pthat = 15 = 882122
  // average xsec = 9.56700e+08
  // no of events in pthat = 30 = 823309
  // average xsec = 7.27938e+07
  // no of events in pthat = 50 = 848247
  // average xsec = 9.84034e+06
  // no of events in pthat = 80 = 834909
  // average xsec = 1.37397e+06
  // no of events in pthat = 120 = 713842
  // average xsec = 2.27780e+05
  // no of events in pthat = 170 = 707705
  // average xsec = 4.28366e+04
  // no of events in pthat = 220 = 707090
  // average xsec = 1.28229e+04
  // no of events in pthat = 280 = 766705
  // average xsec = 3.28043e+03
  // no of events in pthat = 370 = 957400
  // average xsec = 7.57785e+02

  // double pthatWeight[] = {9.56700e+08/882122, 7.27938e+07/823309,
  // 			  9.84034e+06/848247, 1.37397e+06/834909,
  // 			  2.27780e+05/713842, 4.28366e+04/707705,
  // 			  1.28229e+04/707090, 3.28043e+03/766705,
  // 			  7.57785e+02/957400};  
  // return pthatWeight[startfile];
  //! In the new files I made:
  //! pthat 80:
  // number of events: 71557
  // average xsec = 1.32722e+06
  // return (double)1.32722e+06/71557;
  //! In the new files I made:
  //! pthat 170:
  // number of events: 994000
  // average xsec = 
  return (double)44784.9/994000;
  
}

const double pthatbins[] = {15.0, 30.0, 50.0, 80.0, 120.0, 170.0, 220.0, 280.0, 370.0, 99999};

using namespace std;

void makeHistograms_MC(int startfile = 0,
		       int endfile = 1,
		       int radius = 4,
		       std::string algo="",
		       std::string jetType= "PF",
		       std::string kfoutname="pPb_MinBias8TeV_test.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;

  infile_Forest = "EPOSpPb_MinBias_5TeV_forests.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;  

  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 5; 
  TChain * jetTree[N];
  
  string dir[N] = {"hltanalysis", "skimanalysis", Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str()), "hiEvtAnalyzer", "runAnalyzer"};

  string trees[N] = {"HltTree","HltTree","t","HiTree","run"};

  for(int t = 0;t<N;t++){
    jetTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    for(int t = 0;t<N;t++){
      jetTree[t]->Add(filename_Forest.c_str());    
      if(printDebug)cout << "Tree loaded  " << string(dir[t]+"/"+trees[t]).data() << endl;
      if(printDebug)cout << "Entries : " << jetTree[t]->GetEntries() << endl;
    }
  }
  
  for(int i = 0; i<N; ++i)  
    if(i!=2)
      jetTree[2]->AddFriend(jetTree[i]);
  
  int jet40_F;
  int jet60_F;
  int jet80_F;
  int jet40_p_F;
  int jet60_p_F;
  int jet80_p_F;
  int jet40_l1_F;
  int jet60_l1_F;
  int jet80_l1_F;
  int jet40_p_l1_F;
  int jet60_p_l1_F;
  int jet80_p_l1_F;
  // jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v2",&jet40_F);
  // jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v2_Prescl",&jet40_p_F);
  // jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v2",&jet60_F);
  // jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v2_Prescl",&jet60_p_F);
  // jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v2",&jet80_F);
  // jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v2_Prescl",&jet80_p_F);
  // jetTree[0]->SetBranchAddress("L1_SingleJet16_BptxAND_Final",&jet40_l1_F);
  // jetTree[0]->SetBranchAddress("L1_SingleJet16_BptxAND_Prescl",&jet40_p_l1_F);
  // jetTree[0]->SetBranchAddress("L1_SingleJet24_BptxAND_Final",&jet60_l1_F);
  // jetTree[0]->SetBranchAddress("L1_SingleJet24_BptxAND_Prescl",&jet60_p_l1_F);
  // jetTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Final",&jet80_l1_F);
  // jetTree[0]->SetBranchAddress("L1_SingleJet36_BptxAND_Prescl",&jet80_p_l1_F);
  
  int pprimaryvertexFilter_F;
  int pVertexFilterCutGplus_F;
  int pHBHENoiseFilter_F;
  jetTree[1]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);
  jetTree[1]->SetBranchAddress("pPAprimaryVertexFilter",&pprimaryvertexFilter_F);
  jetTree[1]->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus_F);
  
  int nref_F;
  float pt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float m_F[1000];
  float chMax_F[1000];
  int chN_F[1000];
  int neN_F[1000];
  float trkMax_F[1000];
  float chSum_F[1000];
  float phSum_F[1000];
  float neSum_F[1000];
  float trkSum_F[1000];
  float phMax_F[1000];
  float neMax_F[1000];
  float chHardMax_F[1000];
  float trkHardMax_F[1000];
  float chHardSum_F[1000];
  float phHardSum_F[1000];
  float trkHardSum_F[1000];
  float phHardMax_F[1000];
  float eMax_F[1000];
  float muMax_F[1000];
  float eSum_F[1000];
  float muSum_F[1000];
  float jtpu_F[1000];
  float refpt_F[1000];
  float refeta_F[1000];
  int subid_F[1000];
  float refdrjt_F[1000];
  int refparton_F[1000];
  float pthat_F;
  jetTree[2]->SetBranchAddress("nref",&nref_F);
  jetTree[2]->SetBranchAddress("subid",&subid_F);
  jetTree[2]->SetBranchAddress("refdrjt",refdrjt_F);
  jetTree[2]->SetBranchAddress("refparton_flavor",refparton_F);
  jetTree[2]->SetBranchAddress("refpt", &refpt_F);
  jetTree[2]->SetBranchAddress("pthat", &pthat_F);
  jetTree[2]->SetBranchAddress("refeta", &refeta_F);
  jetTree[2]->SetBranchAddress("jtpt",&pt_F);
  jetTree[2]->SetBranchAddress("jteta",&eta_F);
  jetTree[2]->SetBranchAddress("jtphi",&phi_F);
  jetTree[2]->SetBranchAddress("rawpt",&rawpt_F);
  jetTree[2]->SetBranchAddress("jtpu",&jtpu_F);
  jetTree[2]->SetBranchAddress("jtm",&m_F);
  jetTree[2]->SetBranchAddress("chargedMax",&chMax_F);
  jetTree[2]->SetBranchAddress("chargedSum",&chSum_F);
  jetTree[2]->SetBranchAddress("chargedN",&chN_F);
  jetTree[2]->SetBranchAddress("neutralN",&neN_F);
  //jetTree[2]->SetBranchAddress("chargedHardMax",&chMax_F);
  jetTree[2]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetTree[2]->SetBranchAddress("trackMax",&trkMax_F);
  jetTree[2]->SetBranchAddress("trackSum",&trkSum_F);
  //jetTree[2]->SetBranchAddress("trackHardMax",&trkMax_F);
  jetTree[2]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetTree[2]->SetBranchAddress("photonMax",&phMax_F);
  jetTree[2]->SetBranchAddress("photonSum",&phSum_F);
  //jetTree[2]->SetBranchAddress("photonHardMax",&phMax_F);
  jetTree[2]->SetBranchAddress("photonHardSum",&phSum_F);
  jetTree[2]->SetBranchAddress("neutralMax",&neMax_F);
  jetTree[2]->SetBranchAddress("neutralSum",&neSum_F);
  jetTree[2]->SetBranchAddress("eSum",&eSum_F);
  jetTree[2]->SetBranchAddress("eMax",&eMax_F);
  jetTree[2]->SetBranchAddress("muSum",&muSum_F);
  jetTree[2]->SetBranchAddress("muMax",&muMax_F);
  //! pf jet variables
  float jtPfCHF_F[1000];
  float jtPfNHF_F[1000];
  float jtPfCEF_F[1000];
  float jtPfNEF_F[1000];
  float jtPfMUF_F[1000];
  float jtPfCHM_F[1000];
  float jtPfNHM_F[1000];
  float jtPfCEM_F[1000];
  float jtPfNEM_F[1000];
  float jtPfMUM_F[1000];
  jetTree[2]->SetBranchAddress("jtPfCHF",&jtPfCHF_F);
  jetTree[2]->SetBranchAddress("jtPfNHF",&jtPfNHF_F);
  jetTree[2]->SetBranchAddress("jtPfCEF",&jtPfCEF_F);
  jetTree[2]->SetBranchAddress("jtPfNEF",&jtPfNEF_F);
  jetTree[2]->SetBranchAddress("jtPfMUF",&jtPfMUF_F);
  jetTree[2]->SetBranchAddress("jtPfCHM",&jtPfCHM_F);
  jetTree[2]->SetBranchAddress("jtPfNHM",&jtPfNHM_F);
  jetTree[2]->SetBranchAddress("jtPfCEM",&jtPfCEM_F);
  jetTree[2]->SetBranchAddress("jtPfNEM",&jtPfNEM_F);
  jetTree[2]->SetBranchAddress("jtPfMUM",&jtPfMUM_F);

  Float_t vz_F;
  ULong64_t evt_F;
  UInt_t run_F;
  UInt_t lumi_F;
  int hiBin_F;
  int hiNpix_F;
  int hiNtracks_F;
  float hiHF_F;
  float hiHFplus_F;
  float hiHFminus_F;
  float hiZDC_F;
  float hiZDCplus_F;
  float hiZDCminus_F;
  float hiHFhit_F;
  float hiHFhitplus_F;
  float hiHFhitminus_F;
  jetTree[3]->SetBranchAddress("evt",&evt_F);
  jetTree[3]->SetBranchAddress("run",&run_F);
  jetTree[3]->SetBranchAddress("lumi",&lumi_F);
  jetTree[3]->SetBranchAddress("vz",&vz_F);
  jetTree[3]->SetBranchAddress("hiBin",&hiBin_F);
  jetTree[3]->SetBranchAddress("hiHF",&hiHF_F);
  jetTree[3]->SetBranchAddress("hiHFplus",&hiHFplus_F);
  jetTree[3]->SetBranchAddress("hiHFminus",&hiHFminus_F);
  jetTree[3]->SetBranchAddress("hiZDC",&hiZDC_F);
  jetTree[3]->SetBranchAddress("hiZDCplus",&hiZDCplus_F);
  jetTree[3]->SetBranchAddress("hiZDCminus",&hiZDCminus_F);
  jetTree[3]->SetBranchAddress("hiHFhit",&hiHFhit_F);
  jetTree[3]->SetBranchAddress("hiHFhitPlus",&hiHFhitplus_F);
  jetTree[3]->SetBranchAddress("hiHFhitMinus",&hiHFhitminus_F);
  jetTree[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetTree[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);

  Float_t xsec_F;
  jetTree[4]->SetBranchAddress("xsec",&xsec_F);
  
  TFile *fout = new TFile(kfoutname.c_str(),"RECREATE");
  fout->cd();
  
  //! Declare the output histograms to compare data with MC.
  //! Event histograms 
  TH1F * hVz = new TH1F("hVz","",100, -20, 20);
  TH1F * hhiBin = new TH1F("hhiBin","",200, 0, 200);
  TH1F * hNpix = new TH1F("hNpix","", 100, 0, 1200);
  TH1F * hHF = new TH1F("hHF","",500, 0, 500);
  TH1F * hHFplus = new TH1F("hHFplus","",500, 0, 500);
  TH1F * hHFminus = new TH1F("hHFminus","",500, 0, 500);
  TH1F * hZDC = new TH1F("hZDC","",500, 0, 500);
  TH1F * hZDCplus = new TH1F("hZDCplus","",500, 0, 500);
  TH1F * hZDCminus = new TH1F("hZDCminus","",500, 0, 500);
  TH1F * hHFhit = new TH1F("hHFhit","",200, 0, 8000);
  TH1F * hHFhitplus = new TH1F("hHFhitplus","",200, 0, 8000);
  TH1F * hHFhitminus = new TH1F("hHFhitminus","",200, 0, 8000);
  TH1F * hNtracks = new TH1F("hNtracks","", 200, 0, 200);
  TH1F * hpthat = new TH1F("hpthat", "", 500, 0, 1000);
  //! Jet histograms
  TH1F * hjtpt = new TH1F("hjtpt","", 500, 0, 1000);
  TH1F * hrefpt = new TH1F("hrefpt","",500, 0, 1000);
  TH1F * hjteta = new TH1F("hjteta","", 100, -5, 5);
  TH1F * hdijeteta = new TH1F("hdijeteta","",100, -5, 5);
  TH1F * hjtphi = new TH1F("hjtphi","", 100, -4, 4);
  TH1F * hjtm = new TH1F("hjtm","", 50, 0, 50);
  TH1F * hChSumoverRawpT = new TH1F("hChSumoverRawpT","",50, 0, 1.5);
  TH1F * hChMaxoverRawpT = new TH1F("hChMaxoverRawpT","",50, 0, 1.5);
  TH1F * hNeSumoverRawpT = new TH1F("hNeSumoverRawpT","",50, 0, 1.5);
  TH1F * hNeMaxoverRawpT = new TH1F("hNeMaxoverRawpT","",50, 0, 1.5);
  TH1F * hPhSumoverRawpT = new TH1F("hPhSumoverRawpT","",50, 0, 1.5);
  TH1F * hPhMaxoverRawpT = new TH1F("hPhMaxoverRawpT","",50, 0, 1.5);
  TH1F * hElSumoverRawpT = new TH1F("hElSumoverRawpT","",50, 0, 1.5);
  TH1F * hElMaxoverRawpT = new TH1F("hElMaxoverRawpT","",50, 0, 1.5);
  TH1F * hMuSumoverRawpT = new TH1F("hMuSumoverRawpT","",50, 0, 1.5);
  TH1F * hMuMaxoverRawpT = new TH1F("hMuMaxoverRawpT","",50, 0, 1.5);
  TH1F * hCHF = new TH1F("hCHF","", 50, 0, 1.5);
  TH1F * hNHF = new TH1F("hNHF","", 50, 0, 1.5);
  TH1F * hCEF = new TH1F("hCEF","", 50, 0, 1.5);
  TH1F * hNEF = new TH1F("hNEF","", 50, 0, 1.5);
  TH1F * hMUF = new TH1F("hMUF","", 50, 0, 1.5);
  TH1F * hCHM = new TH1F("hCHM","", 70, 0, 70);
  TH1F * hNHM = new TH1F("hNHM","", 70, 0, 70);
  TH1F * hCEM = new TH1F("hCEM","", 70, 0, 70);
  TH1F * hNEM = new TH1F("hNEM","", 70, 0, 70);
  TH1F * hMUM = new TH1F("hMUM","", 70, 0, 70);
  TH1F * hAj = new TH1F("hAj","",50, 0, 1);
  //! 2D histograms with Aj vs different centrality variables
  TH2F * hAjvshiBin = new TH2F("hAjvshiBin","", 200, 0, 200, 50, 0, 1);
  TH2F * hAjvsNpix = new TH2F("hAjvsNpix","", 100, 0, 1200, 50, 0, 1);
  TH2F * hAjvsHF = new TH2F("hAjvsHF","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsHFplus = new TH2F("hAjvsHFplus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsHFminus = new TH2F("hAjvsHFminus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsZDC = new TH2F("hAjvsZDC","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsZDCplus = new TH2F("hAjvsZDCplus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsZDCminus = new TH2F("hAjvsZDCminus","", 500, 0, 500, 50, 0, 1);
  TH2F * hAjvsHFhit = new TH2F("hAjvsHFhit","", 200, 0, 8000, 50, 0, 1);
  TH2F * hAjvsHFhitplus = new TH2F("hAjvsHFhitplus","", 200, 0, 8000, 50, 0, 1);
  TH2F * hAjvsHFhitminus = new TH2F("hAjvsHFhitminus","", 200, 0, 8000, 50, 0, 1);
  TH2F * hAjvsNtracks = new TH2F("hAjvsNtracks","", 200, 0, 200, 50, 0, 1);
  
  //! JES/JER plots 
  TH1F * hJEC[nbins_pt_jec][nbins_eta_jec];
  TH2F * hJEC_Applied[nbins_eta_jec];  
  for(int y = 0; y<nbins_eta_jec; ++y){
    hJEC_Applied[y] = new TH2F(Form("hJEC_Applied_etabin%d", y), "", nbins_pt_jec, ptbins_jec, 60, 0, 3);
    for(int x = 0; x<nbins_pt_jec; ++x){
      hJEC[x][y] = new TH1F(Form("hJEC_ptbin%d_etabin%d",x,y),Form("recopt-genpt/genpt %2.0f < genpt < %2.0f, %2.4f < geneta < %2.4f",ptbins_jec[x], ptbins_jec[x+1], etabins_jec[y], etabins_jec[y+1]),300, 0, 3);
    }
  }
  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetTree[0]->GetEntries();
  if(printDebug) nentries = 100;
  TRandom rnd; 
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    for(int t = 0;t<N;t++)
      jetTree[t]->GetEntry(nEvt);

    if(fabs(vz_F)>15 || !pprimaryvertexFilter_F ||
       !pVertexFilterCutGplus_F || !pHBHENoiseFilter_F)
      continue;
    // if(!jet40_F && !jet60_F && !jet80_F)
    //   continue;
    // if(pthat_F < pthatbins[startfile] || pthat_F >= pthatbins[endfile])
    //   continue;
    
    double ev_weight = 1.0;
    // ev_weight = pthatWeight(startfile);

    double ptlead = 0.0;
    double etalead = 0.0;
    double ptsublead = 0.0;
    double etasublead = 0.0;
    int counter = 0;
    
    for(int jet = 0; jet<nref_F; ++jet){

      if(fabs(eta_F[jet]) > 3.0) continue;
      if(subid_F[jet] != 0) continue;
      // if(pt_F[jet] > 3 * pthat_F) continue;
      if(pt_F[jet] < 20.) continue;
      if(refdrjt_F[jet] > (float)radius/10) continue; 

      double genpt = refpt_F[jet];
      double recpt = pt_F[jet];
      double receta = eta_F[jet];
      double rawpt = rawpt_F[jet];

      if(counter == 0) {
	ptlead = recpt;
	etalead = receta;
      }else if(counter == 1) {
	ptsublead = recpt;
	etasublead = receta;
      }
      counter++;
      
      int etabin = -1;      
      for(int bin = 0; bin<nbins_eta_jec; ++bin){
	if(fabs(receta) > etabins_jec[bin]) etabin = bin;
      }
      if(etabin == -1) continue;
      
      int binx = -1;
      for(int bin = 0; bin<nbins_pt_jec; ++bin){
	if(genpt > ptbins_jec[bin]) binx = bin;
      }
      if(binx == -1) continue;
      
      //! hin JEC:
      hJEC[binx][etabin]->Fill((float)(recpt)/genpt, ev_weight);
      hJEC_Applied[etabin]->Fill(rawpt, (float)recpt/rawpt, ev_weight);      
      hrefpt->Fill(refpt_F[jet], ev_weight);
      hjtpt->Fill(pt_F[jet], ev_weight);
      hjteta->Fill(eta_F[jet], ev_weight);
      hjtphi->Fill(phi_F[jet], ev_weight);
      hjtm->Fill(m_F[jet], ev_weight);
      hChSumoverRawpT->Fill(chSum_F[jet]/rawpt_F[jet], ev_weight);
      hChMaxoverRawpT->Fill(chMax_F[jet]/rawpt_F[jet], ev_weight);
      hNeSumoverRawpT->Fill(neSum_F[jet]/rawpt_F[jet], ev_weight);
      hNeMaxoverRawpT->Fill(neMax_F[jet]/rawpt_F[jet], ev_weight);
      hPhSumoverRawpT->Fill(phSum_F[jet]/rawpt_F[jet], ev_weight);
      hPhMaxoverRawpT->Fill(phMax_F[jet]/rawpt_F[jet], ev_weight);
      hElSumoverRawpT->Fill(eSum_F[jet]/rawpt_F[jet], ev_weight);
      hElMaxoverRawpT->Fill(eMax_F[jet]/rawpt_F[jet], ev_weight);
      hMuSumoverRawpT->Fill(muSum_F[jet]/rawpt_F[jet], ev_weight);
      hMuMaxoverRawpT->Fill(muMax_F[jet]/rawpt_F[jet], ev_weight);
      hCHF->Fill(jtPfCHF_F[jet], ev_weight);
      hNHF->Fill(jtPfNHF_F[jet], ev_weight);
      hCEF->Fill(jtPfCEF_F[jet], ev_weight);
      hNEF->Fill(jtPfNEF_F[jet], ev_weight);
      hMUF->Fill(jtPfMUF_F[jet], ev_weight);
      hCHM->Fill(jtPfCHM_F[jet], ev_weight);
      hNHM->Fill(jtPfNHM_F[jet], ev_weight);
      hCEM->Fill(jtPfCEM_F[jet], ev_weight);
      hNEM->Fill(jtPfNEM_F[jet], ev_weight);
      hMUM->Fill(jtPfMUM_F[jet], ev_weight);
    }

    hpthat->Fill(pthat_F, ev_weight);
    hVz->Fill(vz_F, ev_weight);
    hhiBin->Fill(hiBin_F, ev_weight);
    hNpix->Fill(hiNpix_F, ev_weight);
    hHF->Fill(hiHF_F, ev_weight);
    hHFplus->Fill(hiHFplus_F, ev_weight);
    hHFminus->Fill(hiHFminus_F, ev_weight);
    hZDC->Fill(hiZDC_F, ev_weight);
    hZDCplus->Fill(hiZDCplus_F, ev_weight);
    hZDCminus->Fill(hiZDCminus_F, ev_weight);
    hHFhit->Fill(hiHFhit_F, ev_weight);
    hHFhitplus->Fill(hiHFhitplus_F, ev_weight);
    hHFhitminus->Fill(hiHFhitminus_F, ev_weight);
    hNtracks->Fill(hiNtracks_F, ev_weight);

    double Aj = 0.0;
    double dijeteta = 0.0;
    if(ptlead!=0.0 && ptsublead!=0.0){
      Aj = (double)(ptlead - ptsublead)/(ptlead + ptsublead);
      hAj->Fill(Aj, ev_weight);
      dijeteta = (etalead + etasublead)/2;
      hdijeteta->Fill(dijeteta, ev_weight);
      hAjvshiBin->Fill(hiBin_F, Aj, ev_weight);
      hAjvsNpix->Fill(hiNpix_F, Aj, ev_weight);
      hAjvsNtracks->Fill(hiNtracks_F, Aj, ev_weight);
      hAjvsHF->Fill(hiHF_F, Aj, ev_weight);
      hAjvsHFplus->Fill(hiHFplus_F, Aj, ev_weight);
      hAjvsHFminus->Fill(hiHFminus_F, Aj, ev_weight);
      hAjvsZDC->Fill(hiZDC_F, Aj, ev_weight);
      hAjvsZDCplus->Fill(hiZDCplus_F, Aj, ev_weight);
      hAjvsZDCminus->Fill(hiZDCminus_F, Aj, ev_weight);
      hAjvsHFhit->Fill(hiHFhit_F, Aj, ev_weight);
      hAjvsHFhitplus->Fill(hiHFhitplus_F, Aj, ev_weight);
      hAjvsHFhitminus->Fill(hiHFhitminus_F, Aj, ev_weight);
    }
        
  }// event loop
  
  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
