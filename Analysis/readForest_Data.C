//!
//! Raghav Kunnawalkam Elayavalli
//! Nov 3rd 2016 
//!
//! Macro to read forest files and produce root histograms for the radial moment analysis.
//! Basically we need to plot the first and second radial moment for bachward vs forward
//! and for different pT and Ntrk bins 
//! this takes either MC or Data.
//! for MC there needs to be comparison with gen moments and possbly for unfolding
//! for MC there should also be quark vs gluon 
//! general analysis structure:
//! need to do corrections on the fly! 
//!

#include "readForest.h"

using namespace std;
using namespace fastjet;

int readForest_Data(std::string fileList,
		    int startfile,
		    int endfile,
		    std::string jetAlgo,
		    int radius,
		    std::string jetType,
		    bool debugMode,
		    std::string outfile){
  
  TStopwatch timer;
  timer.Start();
  if(debugMode)std::cout<<std::endl<<"debugMode is ON"<<std::endl;

  std::string dataset = "Data";
  
  std::cout<<"reading filelist "<< fileList << std::endl;
  std::cout<<"Running on " << dataset << std::endl;
  std::cout<<"reading files #'s "<< startfile << " to " << endfile<<std::endl;
  std::cout<<"jet Algo = " << jetAlgo;
  std::cout<<", radius = " << radius;
  std::cout<<", jetType = " << jetType;
  std::cout<<"debugMode = "<<debugMode<<std::endl;
  
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  
  //! declare the histograms 
  TH1F * hJetConst_ungrm[_ntrkbins][_nptbins][_netabins][_nbeta];
  TH1F * hJetConst_grm[_ntrkbins][_nptbins][_netabins][_nbeta];
  TH1F * hJetMass_ungrm[_ntrkbins][_nptbins][_netabins];  
  TH1F * hJetMass_grm[_ntrkbins][_nptbins][_netabins];
  TH2F * hJetMass_vs_pt_ungrm[_ntrkbins][_nptbins][_netabins];
  TH2F * hJetMass_vs_pt_grm[_ntrkbins][_nptbins][_netabins];
  // TH1F * hJetTrack_ungrm[_ntrkbins][_nptbins][_netabins][_nbeta];
  // TH1F * hJetTrack_grm[_ntrkbins][_nptbins][_netabins][_nbeta];
  TH1F * hJetGroomEffect[_ntrkbins][_nptbins][_netabins];    
  TH2F * hJetConstDistribution_ungrm[_ntrkbins][_nptbins][_netabins];
  TH2F * hJetConstDistribution_grm[_ntrkbins][_nptbins][_netabins];
  TH2F * hJetMass_vs_RadMom_ungrm[_ntrkbins][_nptbins][_netabins][_nbeta];
  TH2F * hJetMass_vs_RadMom_grm[_ntrkbins][_nptbins][_netabins][_nbeta];
  TH2F * hJetMoment_vs_pt_ungrm[_ntrkbins][_nptbins][_netabins][_nbeta];
  TH2F * hJetMoment_vs_pt_grm[_ntrkbins][_nptbins][_netabins][_nbeta];  

  std::string trkbins[] = {"0_ntrk_60", "60_ntrk_120", "ntrk_gt120"};
  std::string ptbins[] = {"50_pt_100", "100_pt_300", "pt_gt300"};
  std::string etabins[] = {"m2p6_eta_m0p4", "m0p4_eta_p0p4", "p0p4_eta_p2p6"};
  std::string moments[] = {"firstMoment", "secondMoment"};
  
  for(size_t w = 0; w < _ntrkbins; ++w){
    for(size_t x = 0; x < _nptbins; ++x){
      for(size_t y = 0; y < _netabins; ++y){
	hJetGroomEffect[w][x][y] = new TH1F(Form("hJetGroomEffect_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()),"",100, 0, 1);
	hJetConstDistribution_ungrm[w][x][y] = new TH2F(Form("hJetConstDistribution_ungrm_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 20*radius, 0, radius/10, 100, 0, 1);
	hJetConstDistribution_grm[w][x][y] = new TH2F(Form("hJetConstDistribution_grm_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 20*radius, 0, radius/10, 100, 0, 1);
	hJetMass_ungrm[w][x][y] = new TH1F(Form("hJetMass_ungrm_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 50, 0, 50);
	hJetMass_grm[w][x][y] = new TH1F(Form("hJetMass_grm_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 50, 0, 50);
	hJetMass_vs_pt_ungrm[w][x][y] = new TH2F(Form("hJetMass_vs_pt_ungrm_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 150, 50, 600, 50, 0, 50);
	hJetMass_vs_pt_grm[w][x][y] = new TH2F(Form("hJetMass_vs_pt_grm_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str()), "", 150, 50, 600, 50, 0, 50);
	for(size_t z = 0; z < _nbeta; ++z){
	  hJetConst_ungrm[w][x][y][z] = new TH1F(Form("hJetConst_ungrm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",100, 0, _radMomMax[z]);
	  hJetConst_grm[w][x][y][z] = new TH1F(Form("hJetConst_grm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",100, 0, _radMomMax[z]);
	  // hJetTrack_ungrm[w][x][y][z] = new TH1F(Form("hJetTrack_jetungrm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",100, 0, _radMomMax[z]);
	  // hJetTrack_grm[w][x][y][z] = new TH1F(Form("hJetTrack_jetgrm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",100, 0, _radMomMax[z]);
	  hJetMass_vs_RadMom_ungrm[w][x][y][z] = new TH2F(Form("hJetMass_vs_RadMom_ungrm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",50, 0, 50, 100, 0, _radMomMax[z]);
	  hJetMass_vs_RadMom_grm[w][x][y][z] = new TH2F(Form("hJetMass_vs_RadMom_grm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"",50, 0, 50, 100, 0, _radMomMax[z]);
	  hJetMoment_vs_pt_ungrm[w][x][y][z] = new TH2F(Form("hJetMoment_vs_pt_ungrm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"", 150, 50, 600, 100, 0, _radMomMax[z]);
	  hJetMoment_vs_pt_grm[w][x][y][z] = new TH2F(Form("hJetMoment_vs_pt_grm_%s_%s_%s_%s", trkbins[w].c_str(), ptbins[x].c_str(), etabins[y].c_str(), moments[z].c_str()),"", 150, 50, 600, 100, 0, _radMomMax[z]);
	}
      }
    }
  }
  
  //! setup the analysis framework
  std::ifstream instr_Forest(fileList.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(debugMode)std::cout<<"reading from "<<startfile<<" to "<<endfile<<std::endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }


  TChain * jetTree[_ntrees_Data];

  for(int t = 0;t<_ntrees_Data;t++){
    jetTree[t] = new TChain(std::string(dir_Data[t]).data());
  }//!tree loop ends
  
  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;

    for(int t = 0; t<_ntrees_Data; ++t){
      jetTree[t]->Add(filename_Forest.c_str());    
      if(debugMode)std::cout << "Tree loaded  " << std::string(dir_Data[t]).data() << std::endl;
      if(debugMode)std::cout << "Entries : " << jetTree[t]->GetEntries() << std::endl;
    }
    std::cout<<"Total number of events loaded in HiForest = "<<jetTree[2]->GetEntries()<<std::endl;

  }

  //! friend the trees
  for(int t = 1; t<_ntrees_Data; ++t)
    jetTree[0]->AddFriend(jetTree[t]);

  //! Variables for HiForest
  
  //! jet info
  int nref_F;
  float pt_F[1000];
  float eta_F[1000];
  float phi_F[1000];
  float rawpt_F[1000];
  float jtpu_F[1000];
  float chMax_F[1000];
  float chSum_F[1000];
  int chN_F[1000];
  // float chHardSum_F[1000];
  float trkMax_F[1000];
  float trkSum_F[1000];
  // float trkHardSum_F[1000];
  float phMax_F[1000];
  float phSum_F[1000];
  // float phHardSum_F[1000];
  float neSum_F[1000];
  float neMax_F[1000];
  int neN_F[1000];
  float eMax_F[1000];
  float eSum_F[1000];
  float muMax_F[1000];
  float muSum_F[1000];
  jetTree[0]->SetBranchAddress("nref",&nref_F);
  jetTree[0]->SetBranchAddress("jtpt",pt_F);
  jetTree[0]->SetBranchAddress("jteta",eta_F);
  jetTree[0]->SetBranchAddress("jtphi",phi_F);
  jetTree[0]->SetBranchAddress("rawpt",rawpt_F);
  jetTree[0]->SetBranchAddress("jtpu",jtpu_F);
  jetTree[0]->SetBranchAddress("chargedMax",&chMax_F);
  jetTree[0]->SetBranchAddress("chargedSum",&chSum_F);
  jetTree[0]->SetBranchAddress("chargedN",&chN_F);
  // jetTree[0]->SetBranchAddress("chargedHardSum",&chSum_F);
  jetTree[0]->SetBranchAddress("trackMax",&trkMax_F);
  jetTree[0]->SetBranchAddress("trackSum",&trkSum_F);
  // jetTree[0]->SetBranchAddress("trackHardSum",&trkSum_F);
  jetTree[0]->SetBranchAddress("photonMax",&phMax_F);
  jetTree[0]->SetBranchAddress("photonSum",&phSum_F);
  // jetTree[0]->SetBranchAddress("photonHardSum",&phSum_F);
  jetTree[0]->SetBranchAddress("neutralMax",&neMax_F);
  jetTree[0]->SetBranchAddress("neutralSum",&neSum_F);
  jetTree[0]->SetBranchAddress("neutralN",&neN_F);
  jetTree[0]->SetBranchAddress("eSum",eSum_F);
  jetTree[0]->SetBranchAddress("eMax",eMax_F);
  jetTree[0]->SetBranchAddress("muSum",muSum_F);
  jetTree[0]->SetBranchAddress("muMax",muMax_F);

  //! hlt analysis 
  int minbiasBit, minbiasBit_p;
  int jet40, jet40_p;
  int jet60, jet60_p;
  int jet80, jet80_p;
  int jet100, jet100_p;
  jetTree[1]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_v1", &minbiasBit);
  jetTree[1]->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_v1_Prescl", &minbiasBit_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v3",&jet40);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet40_Eta5p1_v3_Prescl",&jet40_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4",&jet60);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4_Prescl",&jet60_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3",&jet80);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3_Prescl",&jet80_p);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3",&jet100);
  jetTree[1]->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3_Prescl",&jet100_p);
  
  //! noise filters 
  int pBeamScrapingFilter_F;
  int pPAprimaryVertexFilter_F;
  int pHBHENoiseFilter_F;
  jetTree[2]->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter_F);
  jetTree[2]->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter_F);
  jetTree[2]->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&pHBHENoiseFilter_F);

  //! event information 
  ULong64_t evt_F;
  UInt_t run_F;
  UInt_t lumi_F;
  float vz_F;
  int hiNpix_F;
  int hiNpixelTracks_F;
  int hiBin_F;
  float hiHF_F;
  int hiNtracks_F;
  int hiNtracksPtCut_F;
  int hiNtracksEtaCut_F;
  int hiNtracksEtaPtCut_F;
  jetTree[3]->SetBranchAddress("evt",&evt_F);
  jetTree[3]->SetBranchAddress("run",&run_F);
  jetTree[3]->SetBranchAddress("lumi",&lumi_F);
  jetTree[3]->SetBranchAddress("vz",&vz_F);
  jetTree[3]->SetBranchAddress("hiBin",&hiBin_F);
  jetTree[3]->SetBranchAddress("hiHF", &hiHF_F);
  jetTree[3]->SetBranchAddress("hiNpix",&hiNpix_F);
  jetTree[3]->SetBranchAddress("hiNpixelTracks",&hiNpixelTracks_F);
  jetTree[3]->SetBranchAddress("hiNtracks",&hiNtracks_F);
  jetTree[3]->SetBranchAddress("hiNtracksPtCut",&hiNtracksPtCut_F);
  jetTree[3]->SetBranchAddress("hiNtracksEtaCut",&hiNtracksEtaCut_F);
  jetTree[3]->SetBranchAddress("hiNtracksEtaPtCut",&hiNtracksEtaPtCut_F);

  //! Particle flow candidates: vectors
  int nPart_F;
  std::vector<int> *pfID_F = 0;
  std::vector<float> *pfEnergy_F = 0;
  std::vector<float> *pfPt_F = 0;
  std::vector<float> *pfEta_F = 0;
  std::vector<float> *pfPhi_F = 0;
  jetTree[4]->SetBranchAddress("nPFpart", &nPart_F);
  jetTree[4]->SetBranchAddress("pfId", &pfID_F);
  jetTree[4]->SetBranchAddress("pfEnergy", &pfEnergy_F);
  jetTree[4]->SetBranchAddress("pfPt", &pfPt_F);
  jetTree[4]->SetBranchAddress("pfEta", &pfEta_F);
  jetTree[4]->SetBranchAddress("pfPhi", &pfPhi_F);
  
  //! Track trees
  int nTrk_F;
  int trkPt_F[1000];
  int trkEta_F[1000];
  int trkPhi_F[1000];
  //!possibly need to add other high quality cuts
  jetTree[5]->SetBranchAddress("nTrk", &nTrk_F);
  jetTree[5]->SetBranchAddress("trkPt", &trkPt_F);
  jetTree[5]->SetBranchAddress("trkEta", &trkEta_F);
  jetTree[5]->SetBranchAddress("trkPhi", &trkPhi_F);
  

  //! run the analysis: 
  if(debugMode) std::cout<<"Running through all the events now"<<std::endl;
  Long64_t nentries = jetTree[0]->GetEntries();
  if(debugMode) nentries = 1000;
  
  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%10000 == 0)std::cout<<nEvt<<"/"<<nentries<<std::endl;
    if(debugMode){
      std::cout<<"*******************************************"<<std::endl;
      std::cout<<"Start of Event: "<<"  nEvt = "<<nEvt<<std::endl;
    }
    for(int t = 0; t<_ntrees_Data; ++t)
      jetTree[t]->GetEntry(nEvt);

    bool pcollisionEventSelection_F = pBeamScrapingFilter_F && pPAprimaryVertexFilter_F && fabs(vz_F) < 15;
    if(!pcollisionEventSelection_F) continue;

    bool hltSel = jet40 || jet60 || jet80 || jet100;
    if(!hltSel) continue;

    double weight = 1.0;

    //! need to put a track pT cutoff
    int Ntracks = 0;
    for (int  nt = 0; nt < nTrk_F; ++nt) {
      if(trkPt_F[nt]> _trkCut)
	Ntracks++;
    }
    
    int trkbin(-1);
    for (size_t j = 0; j < _ntrkbins; ++j) {
      if (Ntracks >= _trkedges[j] && Ntracks < _trkedges[j+1]) {
	trkbin = j;
      }	
    }
    if(trkbin == -1) continue;

    if(debugMode){
      std::cout<<"Event Multiplicity  = "<<nTrk_F<<"  with chosen tracks no = "<<Ntracks<<"  and bin = "<<trkbin<<std::endl;
      std::cout<<"Event has a total of "<<pfPt_F->size()<<" particle flow candidates"<<std::endl;
      std::cout<<"Jet Multiplicity = "<<nref_F<<", with leading jet pT = "<<pt_F[0]<<std::endl;
    }
    if(pt_F[0] < ptCut) continue;

    //! get the particles from the ParticleFlow collection and make a jet collection from them
    vector<fastjet::PseudoJet> pfObjects;
    for(unsigned ip = 0; ip < pfPt_F->size(); ++ip){
      fastjet::PseudoJet particle(getPX(pfPt_F->at(ip), pfPhi_F->at(ip)),
				  getPY(pfPt_F->at(ip), pfPhi_F->at(ip)),
				  getPZ(pfPt_F->at(ip), pfEta_F->at(ip)),
				  pfEnergy_F->at(ip));
      pfObjects.push_back(particle);
    }

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, radius*0.1);
    fastjet::ClusterSequence cs(pfObjects, jet_def);
    vector<fastjet::PseudoJet> recoJets = fastjet::sorted_by_pt(cs.inclusive_jets(ptCut));

    for (unsigned ij = 0; ij < recoJets.size(); ++ij){
      fastjet::PseudoJet jet = recoJets.at(ij); 

      if(debugMode){
	std::cout<<"Ungroomed Jet "<<std::endl;
	std::cout<<"    raw jet pT = "<<jet.pt()<<", eta = "<<jet.eta()<<", phi = "<<jet.phi()<<std::endl;
      }
      
      //! get the JEC
      vector<JetCorrectorParameters> vpar;   
      FactorizedJetCorrector *JEC = new FactorizedJetCorrector(vpar);
      std::string L2Name;
      std::string L3Name;
      std::string L2L3Res; 
      // L2Name="Summer16_25nsV5_MC_L2Relative_AK4PF.txt";
      // L3Name="Summer16_25nsV5_MC_L3Absolute_AK4PF.txt";
      // L2L3Res="Summer16_25nsV5_MC_L2L3Residual_AK4PF.txt";
      L2Name="Spring16_25nsV8_DATA_L2Relative_AK4PF.txt";
      L3Name="Spring16_25nsV8_DATA_L3Absolute_AK4PF.txt";
      L2L3Res="Spring16_25nsV8_DATA_L2L3Residual_AK4PF.txt";
      // cout << "Using .txt files to update JECs..." << endl;
      // cout << "L2: "<< L2Name << endl;
      // cout << "L3: "<< L3Name << endl;
      JetCorrectorParameters *parl2 = new JetCorrectorParameters(L2Name.c_str());
      JetCorrectorParameters *parl3 = new JetCorrectorParameters(L3Name.c_str());
      JetCorrectorParameters *parl2l3Res = new JetCorrectorParameters(L2L3Res.c_str());
      vpar.push_back(*parl2);
      vpar.push_back(*parl3);
      vpar.push_back(*parl2l3Res);
      JEC = new FactorizedJetCorrector(vpar);
      JEC->setJetEta(jet.eta());
      JEC->setJetPt(jet.pt());
      float jetcorr = JEC->getCorrection();
      float jetpt = jetcorr * jet.pt();

      if(debugMode){
	std::cout<<"   Applying JEC Summer 2016 L2Rel, L3Abs, L2L3Res  correction = "<<jetcorr<<endl;
	std::cout<<"    rec jet pT = "<<jetpt<<", eta = "<<jet.eta()<<", phi = "<<jet.phi()<<std::endl;
      }
      
      int ungrmptbin(-1);
      for (size_t j = 0; j < _nptbins; ++j) {
	if (jetpt >= _ptedges[j] && jetpt < _ptedges[j+1]) {
	  ungrmptbin = j;
	}	
      }
      if(ungrmptbin == -1) continue;

      int etabin(-1);
      for (size_t j = 0; j < _netabins; ++j) {
	if (jet.eta() >= _etaedges[j] && jet.eta() < _etaedges[j+1]) {
	  etabin = j;
	}	
      }
      if(etabin == -1) continue;
	
      hJetMass_ungrm[trkbin][ungrmptbin][etabin]->Fill(jet.m(), weight);

      hJetMass_vs_pt_ungrm[trkbin][ungrmptbin][etabin]->Fill(jetpt, jet.m(), weight);

      for(size_t b = 0; b < _nbeta; ++b){
	double rm = 0.0;
	//! loop over the candidates of the jet and include them in the grid
	for ( unsigned ic = 0; ic<jet.constituents().size(); ++ic ){
	  fastjet::PseudoJet c = jet.constituents().at(ic);
	  double delR = deltaR(jet.eta(), jet.phi_std(), c.eta(), c.phi_std());
	  rm+=(c.pt()*pow(delR,_betaValues[b]))/jetpt;
	  hJetConstDistribution_ungrm[trkbin][ungrmptbin][etabin]->Fill(delR, c.pt()/jetpt, weight);	    
	}
	hJetConst_ungrm[trkbin][ungrmptbin][etabin][b]->Fill(rm, weight);
	hJetMoment_vs_pt_ungrm[trkbin][ungrmptbin][etabin][b]->Fill(jetpt, rm, weight);
	hJetMass_vs_RadMom_ungrm[trkbin][ungrmptbin][etabin][b]->Fill(jet.m(), rm, weight);
	hJetMoment_vs_pt_ungrm[trkbin][ungrmptbin][etabin][b]->Fill(jetpt, rm, weight);
	
	rm = 0.0;
	// foreach (Particle track, tracks.particles()){
	//   double pt = track.momentum().pT();
	//   double delR = deltaR(jet.eta(), jet.phi_std(), track.eta(), track.phi(MINUSPI_PLUSPI));
	//   if (delR <= _jetR){
	//     rm+=(pt*pow(delR,_betaValues[b]))/jet.pt();	      
	//   }
	// }
	// hJetTrack_ungrm[trkbin][ungrmptbin][etabin][b]->Fill(rm, weight);
      }
		
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = fastjet::contrib::RecursiveSymmetryCutBase::scalar_z;
      fastjet::contrib::SoftDrop sd(_beta, _zcut, symmetry_measure, radius*0.1);
      fastjet::PseudoJet sd_jet = sd(jet);    

      if(debugMode){
	std::cout<<"SDGroomed Jet "<<std::endl;
	std::cout<<"    raw jet pT = "<<sd_jet.pt()<<", eta = "<<sd_jet.eta()<<", phi = "<<sd_jet.phi()<<std::endl;
      }

      JEC->setJetEta(sd_jet.eta());
      JEC->setJetPt(sd_jet.pt());
      jetcorr = JEC->getCorrection();
      float sdjetpt = jetcorr * sd_jet.pt();

      if(debugMode){
	std::cout<<"   Applying JEC Summer 2016 L2Rel, L3Abs, L2L3Res  correction = "<<jetcorr<<endl;
	std::cout<<"    rec jet pT = "<<sdjetpt<<", eta = "<<sd_jet.eta()<<", phi = "<<sd_jet.phi()<<std::endl;
      }

      int grmptbin(-1);
      for (size_t j = 0; j < _nptbins; ++j) {
	if (sdjetpt >= _ptedges[j] && sdjetpt < _ptedges[j+1]) {
	  grmptbin = j;
	}	
      }
      if(grmptbin == -1) continue;

      double groomFraction = sd_jet.pt()/jet.pt();
      hJetGroomEffect[trkbin][ungrmptbin][etabin]->Fill(groomFraction, weight);	

      hJetMass_grm[trkbin][ungrmptbin][etabin]->Fill(sd_jet.m(), weight);

      hJetMass_vs_pt_grm[trkbin][ungrmptbin][etabin]->Fill(sdjetpt, jet.m(), weight);
      
      for(size_t b = 0; b < _nbeta; ++b){
	double rm = 0.0;
	//! loop over the candidates of the jet and include them in the grid
	for ( unsigned ic = 0; ic<sd_jet.constituents().size(); ++ic ){
	  fastjet::PseudoJet c = jet.constituents().at(ic);
	  double delR = deltaR(sd_jet.eta(), sd_jet.phi_std(), c.eta(), c.phi_std());
	  rm+=(c.pt()*pow(delR,_betaValues[b]))/sdjetpt;
	  hJetConstDistribution_grm[trkbin][grmptbin][etabin]->Fill(delR, c.pt()/sdjetpt, weight);	    
	}
	hJetConst_grm[trkbin][grmptbin][etabin][b]->Fill(rm, weight);
	hJetMoment_vs_pt_grm[trkbin][grmptbin][etabin][b]->Fill(sdjetpt, rm, weight);
	hJetMass_vs_RadMom_grm[trkbin][grmptbin][etabin][b]->Fill(sd_jet.m(), rm, weight);
	hJetMoment_vs_pt_grm[trkbin][grmptbin][etabin][b]->Fill(sdjetpt, rm, weight);
	rm = 0.0;
	// foreach (Particle track, tracks.particles()){
	//   double pt = track.momentum().pT();
	//   double delR = deltaR(sd_jet.eta(), sd_jet.phi_std(), track.eta(), track.phi(MINUSPI_PLUSPI));
	//   if (delR <= _jetR){
	//     rm+=(pt*pow(delR,_betaValues[b]))/sd_jet.pt();	      
	//   }
	// }
	// hJetTrack_grm[trkbin][grmptbin][etabin][b]->Fill(rm, weight);
      }
    }
    
    //! get the genparticles and make the genjet collections
    // vector<fastjet::PseudoJet> genParticles;    
    
  }//! end of event loop
  
  //! write to output file
  TFile *fout = new TFile(Form("%s", outfile.c_str()),"RECREATE");
  fout->cd();
  for(size_t w = 0; w < _ntrkbins; ++w){
    for(size_t x = 0; x < _nptbins; ++x){
      for(size_t y = 0; y < _netabins; ++y){

	// if(hJetGroomEffect[w][x][y]->GetEntries()>0)
	//   hJetGroomEffect[w][x][y]->Scale(1./hJetGroomEffect[w][x][y]->Integral());
	hJetGroomEffect[w][x][y]->Write();

	// if(hJetConstDistribution_ungrm[w][x][y]->GetEntries()>0)
	//   hJetConstDistribution_ungrm[w][x][y]->Scale(1./hJetConstDistribution_ungrm[w][x][y]->Integral());
	hJetConstDistribution_ungrm[w][x][y]->Write();

	// if(hJetConstDistribution_grm[w][x][y]->GetEntries()>0)
	//   hJetConstDistribution_grm[w][x][y]->Scale(1./hJetConstDistribution_ungrm[w][x][y]->Integral());
	hJetConstDistribution_grm[w][x][y]->Write();
	
	// if(hJetMass_ungrm[w][x][y]->GetEntries()>0)	
	//   hJetMass_ungrm[w][x][y]->Scale(1./hJetMass_ungrm[w][x][y]->Integral());
	hJetMass_ungrm[w][x][y]->Write();

	// if(hJetMass_grm[w][x][y]->GetEntries()>0)	
	//   hJetMass_grm[w][x][y]->Scale(1./hJetMass_grm[w][x][y]->Integral());
	hJetMass_grm[w][x][y]->Write();

	for(size_t z = 0; z < _nbeta; ++z){

	  // if(hJetConst_ungrm[w][x][y][z]->GetEntries()>0)
	  //   hJetConst_ungrm[w][x][y][z]->Scale(1./hJetConst_ungrm[w][x][y][z]->Integral());
	  hJetConst_ungrm[w][x][y][z]->Write();

	  // if(hJetConst_grm[w][x][y][z]->GetEntries()>0)
	  //   hJetConst_grm[w][x][y][z]->Scale(1./hJetConst_grm[w][x][y][z]->Integral());
	  hJetConst_grm[w][x][y][z]->Write();

	  // if(hJetTrack_ungrm[w][x][y][z]->GetEntries()>0)
	  //   hJetTrack_ungrm[w][x][y][z]->Scale(1./hJetTrack_ungrm[w][x][y][z]->Integral());
	  // hJetTrack_ungrm[w][x][y][z]->Write();

	  // if(hJetTrack_grm[w][x][y][z]->GetEntries()>0)
	  //   hJetTrack_grm[w][x][y][z]->Scale(1./hJetTrack_grm[w][x][y][z]->Integral());
	  // hJetTrack_grm[w][x][y][z]->Write();

	  // if(hJetMass_vs_RadMom_ungrm[w][x][y][z]->GetEntries()>0)
	  //   hJetMass_vs_RadMom_ungrm[w][x][y][z]->Scale(1./hJetMass_vs_RadMom_ungrm[w][x][y][z]->Integral());
	  hJetMass_vs_RadMom_ungrm[w][x][y][z]->Write();

	  // if(hJetMass_vs_RadMom_grm[w][x][y][z]->GetEntries()>0)
	  //   hJetMass_vs_RadMom_grm[w][x][y][z]->Scale(1./hJetMass_vs_RadMom_grm[w][x][y][z]->Integral());
	  hJetMass_vs_RadMom_grm[w][x][y][z]->Write();

	  // if(hJetMoment_vs_pt_ungrm[w][x][y][z]->GetEntries()>0)
	  //   hJetMoment_vs_pt_ungrm[w][x][y][z]->Scale(1./hJetMoment_vs_pt_ungrm[w][x][y][z]->Integral());
	  hJetMoment_vs_pt_ungrm[w][x][y][z]->Write();

	  // if(hJetMoment_vs_pt_grm[w][x][y][z]->GetEntries()>0)
	  //   hJetMoment_vs_pt_grm[w][x][y][z]->Scale(1./hJetMoment_vs_pt_grm[w][x][y][z]->Integral());
	  hJetMoment_vs_pt_grm[w][x][y][z]->Write();
	}
      }
    }
  }
  fout->Close();

  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  return 0;
    
}


int main(int argc, char *argv[]){
  
  // error, not enough arguments
  int rStatus = -1;
  if(argc!=9 && argc!=1){
    std::cout<<"for tests on default inputs, do..." <<std::endl;
    std::cout<<"./readForest_Data.exe";
    std::cout<<std::endl<<std::endl;
    std::cout<<"for actually running, do..."<<std::endl;
    std::cout<<"./readForest_Data.exe ";
    std::cout<<"<inputFileList> <startFile> <endFile> ";
    std::cout<<"<jetAlgo> <jetRadius> <jetType> <debugMode> ";
    std::cout<<"<outputFilename> ";
    std::cout<<std::endl<<std::endl;
    std::cout<<"rStatus="<<rStatus<<std::endl;
    return rStatus;
  }
  
  rStatus=1;
  if(argc==1)
    rStatus = readForest_Data();
  else{
    std::string inputFileList=argv[1];
    int startfile= atoi(argv[2]);
    int endfile= atoi(argv[3]);  
    std::string jetAlgo=argv[4];
    int jetRadius= atoi(argv[5]);
    std::string jetType=argv[6];
    bool debug=(bool)atoi(argv[7]);
    std::string outputFileName=argv[8];
    
    rStatus = readForest_Data(inputFileList,
			      startfile,
			      endfile,
			      jetAlgo,
			      jetRadius,
			      jetType,
			      debug,
			      outputFileName);
  }
  std::cout<<"rStatus="<<rStatus<<std::endl;
  return rStatus;
}

