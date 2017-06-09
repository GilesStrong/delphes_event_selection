#include "delphesReader.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"

void delphesReader::Loop() {
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

delphesReader::delphesReader(TTree *tree) : fChain(0) {
	Init(tree);
}

delphesReader::~delphesReader(){
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t delphesReader::GetEntry(Long64_t entry) {
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t delphesReader::LoadTree(Long64_t entry) {
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void delphesReader::Init(TTree *tree) {
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("Event", &Event_, &b_Event_);
	fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
	fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
	fChain->SetBranchAddress("Event.Number", Event_Number, &b_Event_Number);
	fChain->SetBranchAddress("Event.ReadTime", Event_ReadTime, &b_Event_ReadTime);
	fChain->SetBranchAddress("Event.ProcTime", Event_ProcTime, &b_Event_ProcTime);
	fChain->SetBranchAddress("Event.ProcessID", Event_ProcessID, &b_Event_ProcessID);
	fChain->SetBranchAddress("Event.MPI", Event_MPI, &b_Event_MPI);
	fChain->SetBranchAddress("Event.Weight", Event_Weight, &b_Event_Weight);
	fChain->SetBranchAddress("Event.Scale", Event_Scale, &b_Event_Scale);
	fChain->SetBranchAddress("Event.AlphaQED", Event_AlphaQED, &b_Event_AlphaQED);
	fChain->SetBranchAddress("Event.AlphaQCD", Event_AlphaQCD, &b_Event_AlphaQCD);
	fChain->SetBranchAddress("Event.ID1", Event_ID1, &b_Event_ID1);
	fChain->SetBranchAddress("Event.ID2", Event_ID2, &b_Event_ID2);
	fChain->SetBranchAddress("Event.X1", Event_X1, &b_Event_X1);
	fChain->SetBranchAddress("Event.X2", Event_X2, &b_Event_X2);
	fChain->SetBranchAddress("Event.ScalePDF", Event_ScalePDF, &b_Event_ScalePDF);
	fChain->SetBranchAddress("Event.PDF1", Event_PDF1, &b_Event_PDF1);
	fChain->SetBranchAddress("Event.PDF2", Event_PDF2, &b_Event_PDF2);
	fChain->SetBranchAddress("Event_size", &Event_size, &b_Event_size);
	fChain->SetBranchAddress("Particle", &Particle_, &b_Particle_);
	fChain->SetBranchAddress("Particle.fUniqueID", Particle_fUniqueID, &b_Particle_fUniqueID);
	fChain->SetBranchAddress("Particle.fBits", Particle_fBits, &b_Particle_fBits);
	fChain->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
	fChain->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
	fChain->SetBranchAddress("Particle.IsPU", Particle_IsPU, &b_Particle_IsPU);
	fChain->SetBranchAddress("Particle.M1", Particle_M1, &b_Particle_M1);
	fChain->SetBranchAddress("Particle.M2", Particle_M2, &b_Particle_M2);
	fChain->SetBranchAddress("Particle.D1", Particle_D1, &b_Particle_D1);
	fChain->SetBranchAddress("Particle.D2", Particle_D2, &b_Particle_D2);
	fChain->SetBranchAddress("Particle.Charge", Particle_Charge, &b_Particle_Charge);
	fChain->SetBranchAddress("Particle.Mass", Particle_Mass, &b_Particle_Mass);
	fChain->SetBranchAddress("Particle.E", Particle_E, &b_Particle_E);
	fChain->SetBranchAddress("Particle.Px", Particle_Px, &b_Particle_Px);
	fChain->SetBranchAddress("Particle.Py", Particle_Py, &b_Particle_Py);
	fChain->SetBranchAddress("Particle.Pz", Particle_Pz, &b_Particle_Pz);
	fChain->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
	fChain->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
	fChain->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
	fChain->SetBranchAddress("Particle.Rapidity", Particle_Rapidity, &b_Particle_Rapidity);
	fChain->SetBranchAddress("Particle.T", Particle_T, &b_Particle_T);
	fChain->SetBranchAddress("Particle.X", Particle_X, &b_Particle_X);
	fChain->SetBranchAddress("Particle.Y", Particle_Y, &b_Particle_Y);
	fChain->SetBranchAddress("Particle.Z", Particle_Z, &b_Particle_Z);
	fChain->SetBranchAddress("Particle_size", &Particle_size, &b_Particle_size);
	fChain->SetBranchAddress("Track", &Track_, &b_Track_);
	fChain->SetBranchAddress("Track.fUniqueID", Track_fUniqueID, &b_Track_fUniqueID);
	fChain->SetBranchAddress("Track.fBits", Track_fBits, &b_Track_fBits);
	fChain->SetBranchAddress("Track.PID", Track_PID, &b_Track_PID);
	fChain->SetBranchAddress("Track.Charge", Track_Charge, &b_Track_Charge);
	fChain->SetBranchAddress("Track.PT", Track_PT, &b_Track_PT);
	fChain->SetBranchAddress("Track.Eta", Track_Eta, &b_Track_Eta);
	fChain->SetBranchAddress("Track.Phi", Track_Phi, &b_Track_Phi);
	fChain->SetBranchAddress("Track.EtaOuter", Track_EtaOuter, &b_Track_EtaOuter);
	fChain->SetBranchAddress("Track.PhiOuter", Track_PhiOuter, &b_Track_PhiOuter);
	fChain->SetBranchAddress("Track.X", Track_X, &b_Track_X);
	fChain->SetBranchAddress("Track.Y", Track_Y, &b_Track_Y);
	fChain->SetBranchAddress("Track.Z", Track_Z, &b_Track_Z);
	fChain->SetBranchAddress("Track.T", Track_T, &b_Track_T);
	fChain->SetBranchAddress("Track.XOuter", Track_XOuter, &b_Track_XOuter);
	fChain->SetBranchAddress("Track.YOuter", Track_YOuter, &b_Track_YOuter);
	fChain->SetBranchAddress("Track.ZOuter", Track_ZOuter, &b_Track_ZOuter);
	fChain->SetBranchAddress("Track.TOuter", Track_TOuter, &b_Track_TOuter);
	fChain->SetBranchAddress("Track.Dxy", Track_Dxy, &b_Track_Dxy);
	fChain->SetBranchAddress("Track.SDxy", Track_SDxy, &b_Track_SDxy);
	fChain->SetBranchAddress("Track.Xd", Track_Xd, &b_Track_Xd);
	fChain->SetBranchAddress("Track.Yd", Track_Yd, &b_Track_Yd);
	fChain->SetBranchAddress("Track.Zd", Track_Zd, &b_Track_Zd);
	fChain->SetBranchAddress("Track.Particle", Track_Particle, &b_Track_Particle);
	fChain->SetBranchAddress("Track_size", &Track_size, &b_Track_size);
	fChain->SetBranchAddress("Tower", &Tower_, &b_Tower_);
	fChain->SetBranchAddress("Tower.fUniqueID", Tower_fUniqueID, &b_Tower_fUniqueID);
	fChain->SetBranchAddress("Tower.fBits", Tower_fBits, &b_Tower_fBits);
	fChain->SetBranchAddress("Tower.ET", Tower_ET, &b_Tower_ET);
	fChain->SetBranchAddress("Tower.Eta", Tower_Eta, &b_Tower_Eta);
	fChain->SetBranchAddress("Tower.Phi", Tower_Phi, &b_Tower_Phi);
	fChain->SetBranchAddress("Tower.E", Tower_E, &b_Tower_E);
	fChain->SetBranchAddress("Tower.T", Tower_T, &b_Tower_T);
	fChain->SetBranchAddress("Tower.NTimeHits", Tower_NTimeHits, &b_Tower_NTimeHits);
	fChain->SetBranchAddress("Tower.Eem", Tower_Eem, &b_Tower_Eem);
	fChain->SetBranchAddress("Tower.Ehad", Tower_Ehad, &b_Tower_Ehad);
	fChain->SetBranchAddress("Tower.Edges[4]", Tower_Edges, &b_Tower_Edges);
	fChain->SetBranchAddress("Tower.Particles", Tower_Particles, &b_Tower_Particles);
	fChain->SetBranchAddress("Tower_size", &Tower_size, &b_Tower_size);
	fChain->SetBranchAddress("EFlowTrack", &EFlowTrack_, &b_EFlowTrack_);
	fChain->SetBranchAddress("EFlowTrack.fUniqueID", EFlowTrack_fUniqueID, &b_EFlowTrack_fUniqueID);
	fChain->SetBranchAddress("EFlowTrack.fBits", EFlowTrack_fBits, &b_EFlowTrack_fBits);
	fChain->SetBranchAddress("EFlowTrack.PID", EFlowTrack_PID, &b_EFlowTrack_PID);
	fChain->SetBranchAddress("EFlowTrack.Charge", EFlowTrack_Charge, &b_EFlowTrack_Charge);
	fChain->SetBranchAddress("EFlowTrack.PT", EFlowTrack_PT, &b_EFlowTrack_PT);
	fChain->SetBranchAddress("EFlowTrack.Eta", EFlowTrack_Eta, &b_EFlowTrack_Eta);
	fChain->SetBranchAddress("EFlowTrack.Phi", EFlowTrack_Phi, &b_EFlowTrack_Phi);
	fChain->SetBranchAddress("EFlowTrack.EtaOuter", EFlowTrack_EtaOuter, &b_EFlowTrack_EtaOuter);
	fChain->SetBranchAddress("EFlowTrack.PhiOuter", EFlowTrack_PhiOuter, &b_EFlowTrack_PhiOuter);
	fChain->SetBranchAddress("EFlowTrack.X", EFlowTrack_X, &b_EFlowTrack_X);
	fChain->SetBranchAddress("EFlowTrack.Y", EFlowTrack_Y, &b_EFlowTrack_Y);
	fChain->SetBranchAddress("EFlowTrack.Z", EFlowTrack_Z, &b_EFlowTrack_Z);
	fChain->SetBranchAddress("EFlowTrack.T", EFlowTrack_T, &b_EFlowTrack_T);
	fChain->SetBranchAddress("EFlowTrack.XOuter", EFlowTrack_XOuter, &b_EFlowTrack_XOuter);
	fChain->SetBranchAddress("EFlowTrack.YOuter", EFlowTrack_YOuter, &b_EFlowTrack_YOuter);
	fChain->SetBranchAddress("EFlowTrack.ZOuter", EFlowTrack_ZOuter, &b_EFlowTrack_ZOuter);
	fChain->SetBranchAddress("EFlowTrack.TOuter", EFlowTrack_TOuter, &b_EFlowTrack_TOuter);
	fChain->SetBranchAddress("EFlowTrack.Dxy", EFlowTrack_Dxy, &b_EFlowTrack_Dxy);
	fChain->SetBranchAddress("EFlowTrack.SDxy", EFlowTrack_SDxy, &b_EFlowTrack_SDxy);
	fChain->SetBranchAddress("EFlowTrack.Xd", EFlowTrack_Xd, &b_EFlowTrack_Xd);
	fChain->SetBranchAddress("EFlowTrack.Yd", EFlowTrack_Yd, &b_EFlowTrack_Yd);
	fChain->SetBranchAddress("EFlowTrack.Zd", EFlowTrack_Zd, &b_EFlowTrack_Zd);
	fChain->SetBranchAddress("EFlowTrack.Particle", EFlowTrack_Particle, &b_EFlowTrack_Particle);
	fChain->SetBranchAddress("EFlowTrack_size", &EFlowTrack_size, &b_EFlowTrack_size);
	fChain->SetBranchAddress("EFlowPhoton", &EFlowPhoton_, &b_EFlowPhoton_);
	fChain->SetBranchAddress("EFlowPhoton.fUniqueID", EFlowPhoton_fUniqueID, &b_EFlowPhoton_fUniqueID);
	fChain->SetBranchAddress("EFlowPhoton.fBits", EFlowPhoton_fBits, &b_EFlowPhoton_fBits);
	fChain->SetBranchAddress("EFlowPhoton.ET", EFlowPhoton_ET, &b_EFlowPhoton_ET);
	fChain->SetBranchAddress("EFlowPhoton.Eta", EFlowPhoton_Eta, &b_EFlowPhoton_Eta);
	fChain->SetBranchAddress("EFlowPhoton.Phi", EFlowPhoton_Phi, &b_EFlowPhoton_Phi);
	fChain->SetBranchAddress("EFlowPhoton.E", EFlowPhoton_E, &b_EFlowPhoton_E);
	fChain->SetBranchAddress("EFlowPhoton.T", EFlowPhoton_T, &b_EFlowPhoton_T);
	fChain->SetBranchAddress("EFlowPhoton.NTimeHits", EFlowPhoton_NTimeHits, &b_EFlowPhoton_NTimeHits);
	fChain->SetBranchAddress("EFlowPhoton.Eem", EFlowPhoton_Eem, &b_EFlowPhoton_Eem);
	fChain->SetBranchAddress("EFlowPhoton.Ehad", EFlowPhoton_Ehad, &b_EFlowPhoton_Ehad);
	fChain->SetBranchAddress("EFlowPhoton.Edges[4]", EFlowPhoton_Edges, &b_EFlowPhoton_Edges);
	fChain->SetBranchAddress("EFlowPhoton.Particles", EFlowPhoton_Particles, &b_EFlowPhoton_Particles);
	fChain->SetBranchAddress("EFlowPhoton_size", &EFlowPhoton_size, &b_EFlowPhoton_size);
	fChain->SetBranchAddress("EFlowNeutralHadron", &EFlowNeutralHadron_, &b_EFlowNeutralHadron_);
	fChain->SetBranchAddress("EFlowNeutralHadron.fUniqueID", EFlowNeutralHadron_fUniqueID, &b_EFlowNeutralHadron_fUniqueID);
	fChain->SetBranchAddress("EFlowNeutralHadron.fBits", EFlowNeutralHadron_fBits, &b_EFlowNeutralHadron_fBits);
	fChain->SetBranchAddress("EFlowNeutralHadron.ET", EFlowNeutralHadron_ET, &b_EFlowNeutralHadron_ET);
	fChain->SetBranchAddress("EFlowNeutralHadron.Eta", EFlowNeutralHadron_Eta, &b_EFlowNeutralHadron_Eta);
	fChain->SetBranchAddress("EFlowNeutralHadron.Phi", EFlowNeutralHadron_Phi, &b_EFlowNeutralHadron_Phi);
	fChain->SetBranchAddress("EFlowNeutralHadron.E", EFlowNeutralHadron_E, &b_EFlowNeutralHadron_E);
	fChain->SetBranchAddress("EFlowNeutralHadron.T", EFlowNeutralHadron_T, &b_EFlowNeutralHadron_T);
	fChain->SetBranchAddress("EFlowNeutralHadron.NTimeHits", EFlowNeutralHadron_NTimeHits, &b_EFlowNeutralHadron_NTimeHits);
	fChain->SetBranchAddress("EFlowNeutralHadron.Eem", EFlowNeutralHadron_Eem, &b_EFlowNeutralHadron_Eem);
	fChain->SetBranchAddress("EFlowNeutralHadron.Ehad", EFlowNeutralHadron_Ehad, &b_EFlowNeutralHadron_Ehad);
	fChain->SetBranchAddress("EFlowNeutralHadron.Edges[4]", EFlowNeutralHadron_Edges, &b_EFlowNeutralHadron_Edges);
	fChain->SetBranchAddress("EFlowNeutralHadron.Particles", EFlowNeutralHadron_Particles, &b_EFlowNeutralHadron_Particles);
	fChain->SetBranchAddress("EFlowNeutralHadron_size", &EFlowNeutralHadron_size, &b_EFlowNeutralHadron_size);
	fChain->SetBranchAddress("GenJet", &GenJet_, &b_GenJet_);
	fChain->SetBranchAddress("GenJet.fUniqueID", GenJet_fUniqueID, &b_GenJet_fUniqueID);
	fChain->SetBranchAddress("GenJet.fBits", GenJet_fBits, &b_GenJet_fBits);
	fChain->SetBranchAddress("GenJet.PT", GenJet_PT, &b_GenJet_PT);
	fChain->SetBranchAddress("GenJet.Eta", GenJet_Eta, &b_GenJet_Eta);
	fChain->SetBranchAddress("GenJet.Phi", GenJet_Phi, &b_GenJet_Phi);
	fChain->SetBranchAddress("GenJet.T", GenJet_T, &b_GenJet_T);
	fChain->SetBranchAddress("GenJet.Mass", GenJet_Mass, &b_GenJet_Mass);
	fChain->SetBranchAddress("GenJet.DeltaEta", GenJet_DeltaEta, &b_GenJet_DeltaEta);
	fChain->SetBranchAddress("GenJet.DeltaPhi", GenJet_DeltaPhi, &b_GenJet_DeltaPhi);
	fChain->SetBranchAddress("GenJet.Flavor", GenJet_Flavor, &b_GenJet_Flavor);
	fChain->SetBranchAddress("GenJet.FlavorAlgo", GenJet_FlavorAlgo, &b_GenJet_FlavorAlgo);
	fChain->SetBranchAddress("GenJet.FlavorPhys", GenJet_FlavorPhys, &b_GenJet_FlavorPhys);
	fChain->SetBranchAddress("GenJet.BTag", GenJet_BTag, &b_GenJet_BTag);
	fChain->SetBranchAddress("GenJet.BTagAlgo", GenJet_BTagAlgo, &b_GenJet_BTagAlgo);
	fChain->SetBranchAddress("GenJet.BTagPhys", GenJet_BTagPhys, &b_GenJet_BTagPhys);
	fChain->SetBranchAddress("GenJet.TauTag", GenJet_TauTag, &b_GenJet_TauTag);
	fChain->SetBranchAddress("GenJet.Charge", GenJet_Charge, &b_GenJet_Charge);
	fChain->SetBranchAddress("GenJet.EhadOverEem", GenJet_EhadOverEem, &b_GenJet_EhadOverEem);
	fChain->SetBranchAddress("GenJet.NCharged", GenJet_NCharged, &b_GenJet_NCharged);
	fChain->SetBranchAddress("GenJet.NNeutrals", GenJet_NNeutrals, &b_GenJet_NNeutrals);
	fChain->SetBranchAddress("GenJet.Beta", GenJet_Beta, &b_GenJet_Beta);
	fChain->SetBranchAddress("GenJet.BetaStar", GenJet_BetaStar, &b_GenJet_BetaStar);
	fChain->SetBranchAddress("GenJet.MeanSqDeltaR", GenJet_MeanSqDeltaR, &b_GenJet_MeanSqDeltaR);
	fChain->SetBranchAddress("GenJet.PTD", GenJet_PTD, &b_GenJet_PTD);
	fChain->SetBranchAddress("GenJet.FracPt[5]", GenJet_FracPt, &b_GenJet_FracPt);
	fChain->SetBranchAddress("GenJet.Tau[5]", GenJet_Tau, &b_GenJet_Tau);
	fChain->SetBranchAddress("GenJet.TrimmedP4[5]", GenJet_TrimmedP4, &b_GenJet_TrimmedP4);
	fChain->SetBranchAddress("GenJet.PrunedP4[5]", GenJet_PrunedP4, &b_GenJet_PrunedP4);
	fChain->SetBranchAddress("GenJet.SoftDroppedP4[5]", GenJet_SoftDroppedP4, &b_GenJet_SoftDroppedP4);
	fChain->SetBranchAddress("GenJet.NSubJetsTrimmed", GenJet_NSubJetsTrimmed, &b_GenJet_NSubJetsTrimmed);
	fChain->SetBranchAddress("GenJet.NSubJetsPruned", GenJet_NSubJetsPruned, &b_GenJet_NSubJetsPruned);
	fChain->SetBranchAddress("GenJet.NSubJetsSoftDropped", GenJet_NSubJetsSoftDropped, &b_GenJet_NSubJetsSoftDropped);
	fChain->SetBranchAddress("GenJet.Constituents", GenJet_Constituents, &b_GenJet_Constituents);
	fChain->SetBranchAddress("GenJet.Particles", GenJet_Particles, &b_GenJet_Particles);
	fChain->SetBranchAddress("GenJet.Area.fUniqueID", GenJet_Area_fUniqueID, &b_GenJet_Area_fUniqueID);
	fChain->SetBranchAddress("GenJet.Area.fBits", GenJet_Area_fBits, &b_GenJet_Area_fBits);
	fChain->SetBranchAddress("GenJet.Area.fP.fUniqueID", GenJet_Area_fP_fUniqueID, &b_GenJet_Area_fP_fUniqueID);
	fChain->SetBranchAddress("GenJet.Area.fP.fBits", GenJet_Area_fP_fBits, &b_GenJet_Area_fP_fBits);
	fChain->SetBranchAddress("GenJet.Area.fP.fX", GenJet_Area_fP_fX, &b_GenJet_Area_fP_fX);
	fChain->SetBranchAddress("GenJet.Area.fP.fY", GenJet_Area_fP_fY, &b_GenJet_Area_fP_fY);
	fChain->SetBranchAddress("GenJet.Area.fP.fZ", GenJet_Area_fP_fZ, &b_GenJet_Area_fP_fZ);
	fChain->SetBranchAddress("GenJet.Area.fE", GenJet_Area_fE, &b_GenJet_Area_fE);
	fChain->SetBranchAddress("GenJet_size", &GenJet_size, &b_GenJet_size);
	fChain->SetBranchAddress("GenMissingET", &GenMissingET_, &b_GenMissingET_);
	fChain->SetBranchAddress("GenMissingET.fUniqueID", GenMissingET_fUniqueID, &b_GenMissingET_fUniqueID);
	fChain->SetBranchAddress("GenMissingET.fBits", GenMissingET_fBits, &b_GenMissingET_fBits);
	fChain->SetBranchAddress("GenMissingET.MET", GenMissingET_MET, &b_GenMissingET_MET);
	fChain->SetBranchAddress("GenMissingET.Eta", GenMissingET_Eta, &b_GenMissingET_Eta);
	fChain->SetBranchAddress("GenMissingET.Phi", GenMissingET_Phi, &b_GenMissingET_Phi);
	fChain->SetBranchAddress("GenMissingET_size", &GenMissingET_size, &b_GenMissingET_size);
	fChain->SetBranchAddress("Jet", &Jet_, &b_Jet_);
	fChain->SetBranchAddress("Jet.fUniqueID", Jet_fUniqueID, &b_Jet_fUniqueID);
	fChain->SetBranchAddress("Jet.fBits", Jet_fBits, &b_Jet_fBits);
	fChain->SetBranchAddress("Jet.PT", Jet_PT, &b_Jet_PT);
	fChain->SetBranchAddress("Jet.Eta", Jet_Eta, &b_Jet_Eta);
	fChain->SetBranchAddress("Jet.Phi", Jet_Phi, &b_Jet_Phi);
	fChain->SetBranchAddress("Jet.T", Jet_T, &b_Jet_T);
	fChain->SetBranchAddress("Jet.Mass", Jet_Mass, &b_Jet_Mass);
	fChain->SetBranchAddress("Jet.DeltaEta", Jet_DeltaEta, &b_Jet_DeltaEta);
	fChain->SetBranchAddress("Jet.DeltaPhi", Jet_DeltaPhi, &b_Jet_DeltaPhi);
	fChain->SetBranchAddress("Jet.Flavor", Jet_Flavor, &b_Jet_Flavor);
	fChain->SetBranchAddress("Jet.FlavorAlgo", Jet_FlavorAlgo, &b_Jet_FlavorAlgo);
	fChain->SetBranchAddress("Jet.FlavorPhys", Jet_FlavorPhys, &b_Jet_FlavorPhys);
	fChain->SetBranchAddress("Jet.BTag", Jet_BTag, &b_Jet_BTag);
	fChain->SetBranchAddress("Jet.BTagAlgo", Jet_BTagAlgo, &b_Jet_BTagAlgo);
	fChain->SetBranchAddress("Jet.BTagPhys", Jet_BTagPhys, &b_Jet_BTagPhys);
	fChain->SetBranchAddress("Jet.TauTag", Jet_TauTag, &b_Jet_TauTag);
	fChain->SetBranchAddress("Jet.Charge", Jet_Charge, &b_Jet_Charge);
	fChain->SetBranchAddress("Jet.EhadOverEem", Jet_EhadOverEem, &b_Jet_EhadOverEem);
	fChain->SetBranchAddress("Jet.NCharged", Jet_NCharged, &b_Jet_NCharged);
	fChain->SetBranchAddress("Jet.NNeutrals", Jet_NNeutrals, &b_Jet_NNeutrals);
	fChain->SetBranchAddress("Jet.Beta", Jet_Beta, &b_Jet_Beta);
	fChain->SetBranchAddress("Jet.BetaStar", Jet_BetaStar, &b_Jet_BetaStar);
	fChain->SetBranchAddress("Jet.MeanSqDeltaR", Jet_MeanSqDeltaR, &b_Jet_MeanSqDeltaR);
	fChain->SetBranchAddress("Jet.PTD", Jet_PTD, &b_Jet_PTD);
	fChain->SetBranchAddress("Jet.FracPt[5]", Jet_FracPt, &b_Jet_FracPt);
	fChain->SetBranchAddress("Jet.Tau[5]", Jet_Tau, &b_Jet_Tau);
	fChain->SetBranchAddress("Jet.TrimmedP4[5]", Jet_TrimmedP4, &b_Jet_TrimmedP4);
	fChain->SetBranchAddress("Jet.PrunedP4[5]", Jet_PrunedP4, &b_Jet_PrunedP4);
	fChain->SetBranchAddress("Jet.SoftDroppedP4[5]", Jet_SoftDroppedP4, &b_Jet_SoftDroppedP4);
	fChain->SetBranchAddress("Jet.NSubJetsTrimmed", Jet_NSubJetsTrimmed, &b_Jet_NSubJetsTrimmed);
	fChain->SetBranchAddress("Jet.NSubJetsPruned", Jet_NSubJetsPruned, &b_Jet_NSubJetsPruned);
	fChain->SetBranchAddress("Jet.NSubJetsSoftDropped", Jet_NSubJetsSoftDropped, &b_Jet_NSubJetsSoftDropped);
	fChain->SetBranchAddress("Jet.Constituents", Jet_Constituents, &b_Jet_Constituents);
	fChain->SetBranchAddress("Jet.Particles", Jet_Particles, &b_Jet_Particles);
	fChain->SetBranchAddress("Jet.Area.fUniqueID", Jet_Area_fUniqueID, &b_Jet_Area_fUniqueID);
	fChain->SetBranchAddress("Jet.Area.fBits", Jet_Area_fBits, &b_Jet_Area_fBits);
	fChain->SetBranchAddress("Jet.Area.fP.fUniqueID", Jet_Area_fP_fUniqueID, &b_Jet_Area_fP_fUniqueID);
	fChain->SetBranchAddress("Jet.Area.fP.fBits", Jet_Area_fP_fBits, &b_Jet_Area_fP_fBits);
	fChain->SetBranchAddress("Jet.Area.fP.fX", Jet_Area_fP_fX, &b_Jet_Area_fP_fX);
	fChain->SetBranchAddress("Jet.Area.fP.fY", Jet_Area_fP_fY, &b_Jet_Area_fP_fY);
	fChain->SetBranchAddress("Jet.Area.fP.fZ", Jet_Area_fP_fZ, &b_Jet_Area_fP_fZ);
	fChain->SetBranchAddress("Jet.Area.fE", Jet_Area_fE, &b_Jet_Area_fE);
	fChain->SetBranchAddress("Jet_size", &Jet_size, &b_Jet_size);
	fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
	fChain->SetBranchAddress("Electron.fUniqueID", Electron_fUniqueID, &b_Electron_fUniqueID);
	fChain->SetBranchAddress("Electron.fBits", Electron_fBits, &b_Electron_fBits);
	fChain->SetBranchAddress("Electron.PT", Electron_PT, &b_Electron_PT);
	fChain->SetBranchAddress("Electron.Eta", Electron_Eta, &b_Electron_Eta);
	fChain->SetBranchAddress("Electron.Phi", Electron_Phi, &b_Electron_Phi);
	fChain->SetBranchAddress("Electron.T", Electron_T, &b_Electron_T);
	fChain->SetBranchAddress("Electron.Charge", Electron_Charge, &b_Electron_Charge);
	fChain->SetBranchAddress("Electron.EhadOverEem", Electron_EhadOverEem, &b_Electron_EhadOverEem);
	fChain->SetBranchAddress("Electron.Particle", Electron_Particle, &b_Electron_Particle);
	fChain->SetBranchAddress("Electron.IsolationVar", Electron_IsolationVar, &b_Electron_IsolationVar);
	fChain->SetBranchAddress("Electron.IsolationVarRhoCorr", Electron_IsolationVarRhoCorr, &b_Electron_IsolationVarRhoCorr);
	fChain->SetBranchAddress("Electron.SumPtCharged", Electron_SumPtCharged, &b_Electron_SumPtCharged);
	fChain->SetBranchAddress("Electron.SumPtNeutral", Electron_SumPtNeutral, &b_Electron_SumPtNeutral);
	fChain->SetBranchAddress("Electron.SumPtChargedPU", Electron_SumPtChargedPU, &b_Electron_SumPtChargedPU);
	fChain->SetBranchAddress("Electron.SumPt", Electron_SumPt, &b_Electron_SumPt);
	fChain->SetBranchAddress("Electron_size", &Electron_size, &b_Electron_size);
	fChain->SetBranchAddress("Photon", &Photon_, &b_Photon_);
	fChain->SetBranchAddress("Photon.fUniqueID", Photon_fUniqueID, &b_Photon_fUniqueID);
	fChain->SetBranchAddress("Photon.fBits", Photon_fBits, &b_Photon_fBits);
	fChain->SetBranchAddress("Photon.PT", Photon_PT, &b_Photon_PT);
	fChain->SetBranchAddress("Photon.Eta", Photon_Eta, &b_Photon_Eta);
	fChain->SetBranchAddress("Photon.Phi", Photon_Phi, &b_Photon_Phi);
	fChain->SetBranchAddress("Photon.E", Photon_E, &b_Photon_E);
	fChain->SetBranchAddress("Photon.T", Photon_T, &b_Photon_T);
	fChain->SetBranchAddress("Photon.EhadOverEem", Photon_EhadOverEem, &b_Photon_EhadOverEem);
	fChain->SetBranchAddress("Photon.Particles", Photon_Particles, &b_Photon_Particles);
	fChain->SetBranchAddress("Photon.IsolationVar", Photon_IsolationVar, &b_Photon_IsolationVar);
	fChain->SetBranchAddress("Photon.IsolationVarRhoCorr", Photon_IsolationVarRhoCorr, &b_Photon_IsolationVarRhoCorr);
	fChain->SetBranchAddress("Photon.SumPtCharged", Photon_SumPtCharged, &b_Photon_SumPtCharged);
	fChain->SetBranchAddress("Photon.SumPtNeutral", Photon_SumPtNeutral, &b_Photon_SumPtNeutral);
	fChain->SetBranchAddress("Photon.SumPtChargedPU", Photon_SumPtChargedPU, &b_Photon_SumPtChargedPU);
	fChain->SetBranchAddress("Photon.SumPt", Photon_SumPt, &b_Photon_SumPt);
	fChain->SetBranchAddress("Photon_size", &Photon_size, &b_Photon_size);
	fChain->SetBranchAddress("Muon", &Muon_, &b_Muon_);
	fChain->SetBranchAddress("Muon.fUniqueID", Muon_fUniqueID, &b_Muon_fUniqueID);
	fChain->SetBranchAddress("Muon.fBits", Muon_fBits, &b_Muon_fBits);
	fChain->SetBranchAddress("Muon.PT", Muon_PT, &b_Muon_PT);
	fChain->SetBranchAddress("Muon.Eta", Muon_Eta, &b_Muon_Eta);
	fChain->SetBranchAddress("Muon.Phi", Muon_Phi, &b_Muon_Phi);
	fChain->SetBranchAddress("Muon.T", Muon_T, &b_Muon_T);
	fChain->SetBranchAddress("Muon.Charge", Muon_Charge, &b_Muon_Charge);
	fChain->SetBranchAddress("Muon.Particle", Muon_Particle, &b_Muon_Particle);
	fChain->SetBranchAddress("Muon.IsolationVar", Muon_IsolationVar, &b_Muon_IsolationVar);
	fChain->SetBranchAddress("Muon.IsolationVarRhoCorr", Muon_IsolationVarRhoCorr, &b_Muon_IsolationVarRhoCorr);
	fChain->SetBranchAddress("Muon.SumPtCharged", Muon_SumPtCharged, &b_Muon_SumPtCharged);
	fChain->SetBranchAddress("Muon.SumPtNeutral", Muon_SumPtNeutral, &b_Muon_SumPtNeutral);
	fChain->SetBranchAddress("Muon.SumPtChargedPU", Muon_SumPtChargedPU, &b_Muon_SumPtChargedPU);
	fChain->SetBranchAddress("Muon.SumPt", Muon_SumPt, &b_Muon_SumPt);
	fChain->SetBranchAddress("Muon_size", &Muon_size, &b_Muon_size);
	fChain->SetBranchAddress("MissingET", &MissingET_, &b_MissingET_);
	fChain->SetBranchAddress("MissingET.fUniqueID", MissingET_fUniqueID, &b_MissingET_fUniqueID);
	fChain->SetBranchAddress("MissingET.fBits", MissingET_fBits, &b_MissingET_fBits);
	fChain->SetBranchAddress("MissingET.MET", MissingET_MET, &b_MissingET_MET);
	fChain->SetBranchAddress("MissingET.Eta", MissingET_Eta, &b_MissingET_Eta);
	fChain->SetBranchAddress("MissingET.Phi", MissingET_Phi, &b_MissingET_Phi);
	fChain->SetBranchAddress("MissingET_size", &MissingET_size, &b_MissingET_size);
	fChain->SetBranchAddress("ScalarHT", &ScalarHT_, &b_ScalarHT_);
	fChain->SetBranchAddress("ScalarHT.fUniqueID", ScalarHT_fUniqueID, &b_ScalarHT_fUniqueID);
	fChain->SetBranchAddress("ScalarHT.fBits", ScalarHT_fBits, &b_ScalarHT_fBits);
	fChain->SetBranchAddress("ScalarHT.HT", ScalarHT_HT, &b_ScalarHT_HT);
	fChain->SetBranchAddress("ScalarHT_size", &ScalarHT_size, &b_ScalarHT_size);
	Notify();
}

Bool_t delphesReader::Notify() {
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void delphesReader::Show(Long64_t entry) {
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t delphesReader::Cut(Long64_t entry) {
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}
