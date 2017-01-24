//! Raghav Kunnawalkam Elayavalli
//! Jan 12th 2017
//! Rutgers 
//! for questions or comments: raghav.k.e at CERN dot CH

//! Macro to read Data hiforest and make simple histograms 

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

double eventPrescl(int *trgDec, int *treePrescl, double triggerPt){
  float weight = 0.0;
  if(trgDec[0] && triggerPt>=40.  && triggerPt<60. )
    weight = treePrescl[0];
  if(trgDec[1] && triggerPt>=60.  && triggerPt<80. )
    weight = treePrescl[1];
  if(trgDec[2] && triggerPt>=80. )
    weight = treePrescl[2];
  return weight;
}

int _centedges[] = {0, 40, 80, 120, 200};
const int _ncentbins = sizeof(_centedges)/sizeof(int)-1;
std::string centbins[] = {"0_cent_20", "20_cent_40", "40_cent_60", "60_cent_100"};

double _ptedges[] = {50., 100., 300., 500, 2000.};
const int _nptbins = sizeof(_ptedges)/sizeof(double)-1;
std::string ptbins[] = {"50_pt_100", "100_pt_300", "300_pt_500", "pt_gt500"};
  
double _etaedges[] = {-3.5, -2.0, -0.5, 0.5, 2.0, 3.5};
const int _netabins = sizeof(_etaedges)/sizeof(double)-1;
std::string etabins[] = {"m3p5_eta_m2p0", "m2p0_eta_m0p5",
			 "m0p5_eta_p0p5",
			 "p0p5_eta_p2p0", "p2p0_eta_3p5"};

double _betaValues[] = {0.5, 1., 2., 3.};
double _radMomMax[] = {0.8, 0.5, 0.3, 0.05};
const int _nbeta = sizeof(_betaValues)/sizeof(double);
std::string moments[] = {"halfMoment", "firstMoment", "secondMoment", "thirdMoment"};
    
using namespace std;

void makeHistograms_Data(int startfile = 0,
			 int endfile = 10,
			 int radius = 4,
			 std::string algo="",
			 std::string jetType= "PF",
			 std::string kfoutname="pPb_5p02TeV_MinBias_ak4PF.root"){
  
  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  if(printDebug)cout<<"radius = "<<radius<<endl;
  
  TDatime date;

  std::string infile_Forest;
  infile_Forest = "pPb_MinBias5TeV_forests.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 7; 

  TChain * jetTree[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = Form("ak%s%d%sJetAnalyzer", algo.c_str(), radius, jetType.c_str());
  dir[3] = "hiEvtAnalyzer";
  dir[4] = "hltobject";
  dir[5] = "hltobject";
  dir[6] = "hltobject";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "t",
    "HiTree",
    "HLT_PAAK4PFJet40_Eta5p1_v",
    "HLT_PAAK4PFJet60_Eta5p1_v",
    "HLT_PAAK4PFJet80_Eta5p1_v"
  };

  for(int t = 0;t<N;t++){
    jetTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends

  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    for(int i = 0; i<N; ++i){
      jetTree[i]->Add(filename_Forest.c_str());    
      if(printDebug)cout << "Tree loaded  " << string(dir[i]+"/"+trees[i]).data() << endl;
      if(printDebug)cout << "Entries : " << jetTree[i]->GetEntries() << endl;
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
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v3",&jet40_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v3_Prescl",&jet40_p_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4",&jet60_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4_Prescl",&jet60_p_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3",&jet80_F);
  jetTree[0]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3_Prescl",&jet80_p_F);
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
  float jtpu_F[1000];
  float jtarea_F[1000];
  int jtnCands_F[1000];
  int jtnChCands_F[1000];
  int jtnNeCands_F[1000];
  float jtMByPt_F[1000];
  float jtRMSCand_F[1000];
  float jtAxis1_F[1000];
  float jtAxis2_F[1000];
  float jtSigma_F[1000];
  float jtR_F[1000];
  float jtpTD_F[1000];
  float jtrm0p5_F[1000];
  float jtrm1_F[1000];
  float jtrm2_F[1000];
  float jtrm3_F[1000];
  float jtpull_F[1000];
  float jtSDm_F[1000];
  float jtSDpt_F[1000];
  float jtSDptFrac_F[1000];
  float jtSDrm0p5_F[1000];
  float jtSDrm1_F[1000];
  float jtSDrm2_F[1000];
  float jtSDrm3_F[1000];
  jetTree[2]->SetBranchAddress("nref",&nref_F);
  jetTree[2]->SetBranchAddress("jtpt",&pt_F);
  jetTree[2]->SetBranchAddress("jteta",&eta_F);
  jetTree[2]->SetBranchAddress("jtphi",&phi_F);
  jetTree[2]->SetBranchAddress("rawpt",&rawpt_F);
  jetTree[2]->SetBranchAddress("jtpu",&jtpu_F);
  jetTree[2]->SetBranchAddress("jtm",&m_F);
  jetTree[2]->SetBranchAddress("jtarea",&jtarea_F);
  /*
  jetTree[2]->SetBranchAddress("jtnCands",&jtnCands_F);
  jetTree[2]->SetBranchAddress("jtnChCands",&jtnChCands_F);
  jetTree[2]->SetBranchAddress("jtnNeCands",&jtnNeCands_F);
  jetTree[2]->SetBranchAddress("jtMByPt",&jtMByPt_F);
  jetTree[2]->SetBranchAddress("jtRMSCand",&jtRMSCand_F);
  jetTree[2]->SetBranchAddress("jtAxis1",&jtAxis1_F);
  jetTree[2]->SetBranchAddress("jtAxis2",&jtAxis2_F);
  jetTree[2]->SetBranchAddress("jtSigma",&jtSigma_F);
  jetTree[2]->SetBranchAddress("jtR",&jtR_F);
  jetTree[2]->SetBranchAddress("jtpTD",&jtpTD_F);
  jetTree[2]->SetBranchAddress("jtpull",&jtpull_F);
  jetTree[2]->SetBranchAddress("jtrm0p5",&jtrm0p5_F);
  jetTree[2]->SetBranchAddress("jtrm1",&jtrm1_F);
  jetTree[2]->SetBranchAddress("jtrm2",&jtrm2_F);
  jetTree[2]->SetBranchAddress("jtrm3",&jtrm3_F);
  jetTree[2]->SetBranchAddress("jtSDm",&jtSDm_F);
  jetTree[2]->SetBranchAddress("jtSDpt",&jtSDpt_F);
  jetTree[2]->SetBranchAddress("jtSDptFrac",&jtSDptFrac_F);
  jetTree[2]->SetBranchAddress("jtSDrm0p5",&jtSDrm0p5_F);
  jetTree[2]->SetBranchAddress("jtSDrm1",&jtSDrm1_F);
  jetTree[2]->SetBranchAddress("jtSDrm2",&jtSDrm2_F);
  jetTree[2]->SetBranchAddress("jtSDrm3",&jtSDrm3_F);
  */
  
  //! pf jet variables
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

  //! vectors for hltobject tree
  std::vector<double>  *trgObjpt_40, *trgObjpt_60, *trgObjpt_80;
  jetTree[4]->SetBranchAddress("pt",&trgObjpt_40);
  jetTree[5]->SetBranchAddress("pt",&trgObjpt_60);  
  jetTree[6]->SetBranchAddress("pt",&trgObjpt_80);  
  
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
  //! Jet histograms
  TH1F * hjtpt = new TH1F("hjtpt","", 500, 0, 1000);
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

  /*
  //! declare the histograms 
  TH1F * hJetConst_ungrm[_ncentbins][_nptbins][_netabins][_nbeta];  
  TH1F * hJetConst_grm[_ncentbins][_nptbins][_netabins][_nbeta];
  TH1F * hJetMass_ungrm[_ncentbins][_nptbins][_netabins];  
  TH1F * hJetMass_grm[_ncentbins][_nptbins][_netabins];
  TH2F * hJetMass_vs_pt_ungrm[_ncentbins][_nptbins][_netabins];
  TH2F * hJetMass_vs_pt_grm[_ncentbins][_nptbins][_netabins];  
  TH1F * hJetGroomEffect[_ncentbins][_nptbins][_netabins];    
  // TH2F * hJetConstDistribution_ungrm[_ncentbins][_nptbins][_netabins];
  // TH2F * hJetConstDistribution_grm[_ncentbins][_nptbins][_netabins];
  TH2F * hJetMass_vs_RadMom_ungrm[_ncentbins][_nptbins][_netabins][_nbeta];
  TH2F * hJetMass_vs_RadMom_grm[_ncentbins][_nptbins][_netabins][_nbeta];
  TH2F * hJetMoment_vs_pt_ungrm[_ncentbins][_nptbins][_netabins][_nbeta];
  TH2F * hJetMoment_vs_pt_grm[_ncentbins][_nptbins][_netabins][_nbeta];  

  for(size_t w = 0; w < _ncentbins; ++w){
    for(size_t x = 0; x < _nptbins; ++x){
      for(size_t y = 0; y < _netabins; ++y){
	hJetGroomEffect[w][x][y] = new TH1F(Form("hJetGroomEffect_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()),"",100, 0, 1);
	// hJetConstDistribution_ungrm[w][x][y] = new TH2F(Form("hJetConstDistribution_ungrm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 20*radius, 0, radius/10, 100, 0, 1);
	// hJetConstDistribution_grm[w][x][y] = new TH2F(Form("hJetConstDistribution_grm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 20*radius, 0, radius/10, 100, 0, 1);
	hJetMass_ungrm[w][x][y] = new TH1F(Form("hJetMass_ungrm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 50, 0, 50);
	hJetMass_grm[w][x][y] = new TH1F(Form("hJetMass_grm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 50, 0, 50);
	hJetMass_vs_pt_ungrm[w][x][y] = new TH2F(Form("hJetMass_vs_pt_ungrm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 150, 50, 600, 50, 0, 50);
	hJetMass_vs_pt_grm[w][x][y] = new TH2F(Form("hJetMass_vs_pt_grm_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 150, 50, 600, 50, 0, 50);
	for(size_t z = 0; z < _nbeta; ++z){
	  hJetConst_ungrm[w][x][y][z] = new TH1F(Form("hJetConst_ungrm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",100, 0, _radMomMax[z]);
	  hJetConst_grm[w][x][y][z] = new TH1F(Form("hJetConst_grm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",100, 0, _radMomMax[z]);
	  hJetMass_vs_RadMom_ungrm[w][x][y][z] = new TH2F(Form("hJetMass_vs_RadMom_ungrm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"", 100, 0, _radMomMax[z],50, 0, 50);
	  hJetMass_vs_RadMom_grm[w][x][y][z] = new TH2F(Form("hJetMass_vs_RadMom_grm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"", 100, 0, _radMomMax[z],50, 0, 50);
	  hJetMoment_vs_pt_ungrm[w][x][y][z] = new TH2F(Form("hJetMoment_vs_pt_ungrm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"", 150, 50, 600, 100, 0, _radMomMax[z]);
	  hJetMoment_vs_pt_grm[w][x][y][z] = new TH2F(Form("hJetMoment_vs_pt_grm_%s_%s_%s_%s", centbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"", 150, 50, 600, 100, 0, _radMomMax[z]);
	}
      }
    }
  }
  */
  
  //! now start the event loop for each file.  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetTree[0]->GetEntries();
  if(printDebug) nentries = 100;
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;

    for(int t = 0;t<N;t++)
      jetTree[t]->GetEntry(nEvt);
    
    //! primary event cleaning
    if(fabs(vz_F)>15 || !pprimaryvertexFilter_F ||
       !pVertexFilterCutGplus_F || !pHBHENoiseFilter_F)
      continue;
    if(pt_F[0]<20) continue;
    // if(!jet40_F && !jet60_F && !jet80_F)
    //   continue;
    
    double ev_weight = 1.0;
    // double maxTrigpT = 0.0;
    // if(trgObjpt_40->size()>0){
    //   for(unsigned it = 0; it<trgObjpt_40->size(); ++it){
    // 	if(maxTrigpT < trgObjpt_40->at(it))
    // 	  maxTrigpT = trgObjpt_40->at(it);
    //   }
    // }
    // if(trgObjpt_60->size()>0){
    //   for(unsigned it = 0; it<trgObjpt_60->size(); ++it){
    // 	if(maxTrigpT < trgObjpt_60->at(it))
    // 	  maxTrigpT = trgObjpt_60->at(it);
    //   }
    // }
    // if(trgObjpt_80->size()>0){
    //   for(unsigned it = 0; it<trgObjpt_80->size(); ++it){
    // 	if(maxTrigpT < trgObjpt_80->at(it))
    // 	  maxTrigpT = trgObjpt_80->at(it);
    //   }
    // }

    // int presclDes[] = {jet40_F, jet60_F, jet80_F};
    // int presclWgt[] = {jet40_p_F, jet60_p_F, jet80_p_F};
    
    // ev_weight = eventPrescl(presclDes, presclWgt, maxTrigpT);
    // if(printDebug){
    //   cout<<"Prescl weight = "<<ev_weight<<endl;
    // }
    // ev_weight = 1.0;

    //! fill in the event histograms
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

    double ptlead = 0.0;
    double etalead = 0.0;
    double ptsublead = 0.0;
    double etasublead = 0.0;
    int counter = 0;

    int centbin = -1;
    for(int bin = 0; bin<_ncentbins; ++bin){
      if(hiBin_F >= _centedges[bin] && hiBin_F < _centedges[bin+1])
	centbin = bin;
    }
    if(centbin == -1)
      continue;
    
    //! Loop over the jets 
    for(int jet = 0; jet<nref_F; ++jet){
      //! no JetID at the moment
      if(fabs(eta_F[jet])>3.0) continue;
      if(pt_F[jet] < 20.) continue;

      if(counter == 0) {
	ptlead = pt_F[jet];
	etalead = eta_F[jet];
      }else if(counter == 1) {
	ptsublead = pt_F[jet];
	etasublead = eta_F[jet];
      }
      counter++;
      
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

      int ptbin = -1;
      for(int bin = 0; bin<_nptbins; ++bin){
	if(pt_F[jet] >= _ptedges[bin] && pt_F[jet] < _ptedges[bin+1])
	  ptbin = bin;
      }
      if(ptbin == -1)
	continue;

      int etabin = -1;
      for(int bin = 0; bin<_netabins; ++bin){
	if(eta_F[jet] >= _etaedges[bin] && eta_F[jet] < _etaedges[bin+1])
	  etabin = bin;
      }
      if(etabin == -1)
	continue;      

      /*
      hJetMass_ungrm[centbin][ptbin][etabin]->Fill(m_F[jet], ev_weight);  
      hJetMass_grm[centbin][ptbin][etabin]->Fill(jtSDm_F[jet], ev_weight);
      hJetMass_vs_pt_ungrm[centbin][ptbin][etabin]->Fill(pt_F[jet], m_F[jet], ev_weight);
      hJetMass_vs_pt_grm[centbin][ptbin][etabin]->Fill(pt_F[jet], jtSDm_F[jet], ev_weight);  
      hJetGroomEffect[centbin][ptbin][etabin]->Fill(jtSDpt_F[jet]/pt_F[jet], ev_weight);    
      // hJetConstDistribution_ungrm[centbin][ptbin][etabin]->Fill(, ev_weight);
      // hJetConstDistribution_grm[centbin][ptbin][etabin]->Fill(, ev_weight);

      //! 0.5 radial moment
      hJetConst_ungrm[centbin][ptbin][etabin][0]->Fill(jtrm0p5_F[jet], ev_weight);  
      hJetConst_grm[centbin][ptbin][etabin][0]->Fill(jtSDrm0p5_F[jet], ev_weight);
      hJetMass_vs_RadMom_ungrm[centbin][ptbin][etabin][0]->Fill(jtrm0p5_F[jet], m_F[jet], ev_weight);
      hJetMass_vs_RadMom_grm[centbin][ptbin][etabin][0]->Fill(jtSDrm0p5_F[jet], jtSDm_F[jet], ev_weight);
      hJetMoment_vs_pt_ungrm[centbin][ptbin][etabin][0]->Fill(pt_F[jet], jtrm0p5_F[jet], ev_weight);
      hJetMoment_vs_pt_grm[centbin][ptbin][etabin][0]->Fill(pt_F[jet], jtSDrm0p5_F[jet], ev_weight);        
      //! 1 radial moment
      hJetConst_ungrm[centbin][ptbin][etabin][1]->Fill(jtrm1_F[jet], ev_weight);  
      hJetConst_grm[centbin][ptbin][etabin][1]->Fill(jtSDrm1_F[jet], ev_weight);
      hJetMass_vs_RadMom_ungrm[centbin][ptbin][etabin][1]->Fill(jtrm1_F[jet], m_F[jet], ev_weight);
      hJetMass_vs_RadMom_grm[centbin][ptbin][etabin][1]->Fill(jtSDrm1_F[jet], jtSDm_F[jet], ev_weight);
      hJetMoment_vs_pt_ungrm[centbin][ptbin][etabin][1]->Fill(pt_F[jet], jtrm1_F[jet], ev_weight);
      hJetMoment_vs_pt_grm[centbin][ptbin][etabin][1]->Fill(pt_F[jet], jtSDrm1_F[jet], ev_weight);        
      //! 2 radial moment
      hJetConst_ungrm[centbin][ptbin][etabin][2]->Fill(jtrm2_F[jet], ev_weight);  
      hJetConst_grm[centbin][ptbin][etabin][2]->Fill(jtSDrm2_F[jet], ev_weight);
      hJetMass_vs_RadMom_ungrm[centbin][ptbin][etabin][2]->Fill(jtrm2_F[jet], m_F[jet], ev_weight);
      hJetMass_vs_RadMom_grm[centbin][ptbin][etabin][2]->Fill(jtSDrm2_F[jet], jtSDm_F[jet], ev_weight);
      hJetMoment_vs_pt_ungrm[centbin][ptbin][etabin][2]->Fill(pt_F[jet], jtrm2_F[jet], ev_weight);
      hJetMoment_vs_pt_grm[centbin][ptbin][etabin][2]->Fill(pt_F[jet], jtSDrm2_F[jet], ev_weight);        
      //! 3 radial moment
      hJetConst_ungrm[centbin][ptbin][etabin][3]->Fill(jtrm3_F[jet], ev_weight);  
      hJetConst_grm[centbin][ptbin][etabin][3]->Fill(jtSDrm3_F[jet], ev_weight);
      hJetMass_vs_RadMom_ungrm[centbin][ptbin][etabin][3]->Fill(jtrm3_F[jet], m_F[jet], ev_weight);
      hJetMass_vs_RadMom_grm[centbin][ptbin][etabin][3]->Fill(jtSDrm3_F[jet], jtSDm_F[jet], ev_weight);
      hJetMoment_vs_pt_ungrm[centbin][ptbin][etabin][3]->Fill(pt_F[jet], jtrm3_F[jet], ev_weight);
      hJetMoment_vs_pt_grm[centbin][ptbin][etabin][3]->Fill(pt_F[jet], jtSDrm3_F[jet], ev_weight);        
      */
    }// jet loop
    
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
    
    if(printDebug)cout<<endl;

  }// event loop

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
}//macro end
