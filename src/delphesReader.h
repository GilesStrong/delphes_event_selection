//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 20 14:32:55 2018 by ROOT version 6.06/06
// from TTree Delphes/Analysis tree
// found on file: GluGluToHHTo2B2Tau_node_2_14TeV-madgraph_1_0.root
//////////////////////////////////////////////////////////

#ifndef delphesReader_h
#define delphesReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TObject.h>
#include "TClonesArray.h"
#include "TObject.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxEvent = 1;
const Int_t kMaxWeight = 919;
const Int_t kMaxParticle = 611;
const Int_t kMaxVertex = 258;
const Int_t kMaxGenJet = 17;
const Int_t kMaxGenJetAK8 = 5;
const Int_t kMaxGenMissingET = 1;
const Int_t kMaxPhotonLoose = 7;
const Int_t kMaxPhotonTight = 5;
const Int_t kMaxElectron = 4;
const Int_t kMaxMuonLoose = 14;
const Int_t kMaxMuonTight = 4;
const Int_t kMaxElectronCHS = 29;
const Int_t kMaxMuonLooseCHS = 14;
const Int_t kMaxMuonTightCHS = 5;
const Int_t kMaxJet = 28;
const Int_t kMaxJetPUPPI = 15;
const Int_t kMaxJetAK8 = 15;
const Int_t kMaxJetPUPPIAK8 = 5;
const Int_t kMaxRho = 5;
const Int_t kMaxMissingET = 1;
const Int_t kMaxPuppiMissingET = 1;
const Int_t kMaxGenPileUpMissingET = 1;
const Int_t kMaxScalarHT = 1;

class delphesReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   Long64_t        Event_Number[kMaxEvent];   //[Event_]
   Float_t         Event_ReadTime[kMaxEvent];   //[Event_]
   Float_t         Event_ProcTime[kMaxEvent];   //[Event_]
   Int_t           Event_ProcessID[kMaxEvent];   //[Event_]
   Int_t           Event_MPI[kMaxEvent];   //[Event_]
   Float_t         Event_Weight[kMaxEvent];   //[Event_]
   Float_t         Event_Scale[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQED[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQCD[kMaxEvent];   //[Event_]
   Int_t           Event_ID1[kMaxEvent];   //[Event_]
   Int_t           Event_ID2[kMaxEvent];   //[Event_]
   Float_t         Event_X1[kMaxEvent];   //[Event_]
   Float_t         Event_X2[kMaxEvent];   //[Event_]
   Float_t         Event_ScalePDF[kMaxEvent];   //[Event_]
   Float_t         Event_PDF1[kMaxEvent];   //[Event_]
   Float_t         Event_PDF2[kMaxEvent];   //[Event_]
   Int_t           Event_size;
   Int_t           Weight_;
   UInt_t          Weight_fUniqueID[kMaxWeight];   //[Weight_]
   UInt_t          Weight_fBits[kMaxWeight];   //[Weight_]
   Float_t         Weight_Weight[kMaxWeight];   //[Weight_]
   Int_t           Weight_size;
   Int_t           Particle_;
   UInt_t          Particle_fUniqueID[kMaxParticle];   //[Particle_]
   UInt_t          Particle_fBits[kMaxParticle];   //[Particle_]
   Int_t           Particle_PID[kMaxParticle];   //[Particle_]
   Int_t           Particle_Status[kMaxParticle];   //[Particle_]
   Int_t           Particle_IsPU[kMaxParticle];   //[Particle_]
   Int_t           Particle_M1[kMaxParticle];   //[Particle_]
   Int_t           Particle_M2[kMaxParticle];   //[Particle_]
   Int_t           Particle_D1[kMaxParticle];   //[Particle_]
   Int_t           Particle_D2[kMaxParticle];   //[Particle_]
   Int_t           Particle_Charge[kMaxParticle];   //[Particle_]
   Float_t         Particle_Mass[kMaxParticle];   //[Particle_]
   Float_t         Particle_E[kMaxParticle];   //[Particle_]
   Float_t         Particle_Px[kMaxParticle];   //[Particle_]
   Float_t         Particle_Py[kMaxParticle];   //[Particle_]
   Float_t         Particle_Pz[kMaxParticle];   //[Particle_]
   Float_t         Particle_P[kMaxParticle];   //[Particle_]
   Float_t         Particle_PT[kMaxParticle];   //[Particle_]
   Float_t         Particle_Eta[kMaxParticle];   //[Particle_]
   Float_t         Particle_Phi[kMaxParticle];   //[Particle_]
   Float_t         Particle_Rapidity[kMaxParticle];   //[Particle_]
   Float_t         Particle_CtgTheta[kMaxParticle];   //[Particle_]
   Float_t         Particle_D0[kMaxParticle];   //[Particle_]
   Float_t         Particle_DZ[kMaxParticle];   //[Particle_]
   Float_t         Particle_T[kMaxParticle];   //[Particle_]
   Float_t         Particle_X[kMaxParticle];   //[Particle_]
   Float_t         Particle_Y[kMaxParticle];   //[Particle_]
   Float_t         Particle_Z[kMaxParticle];   //[Particle_]
   Int_t           Particle_size;
   Int_t           Vertex_;
   UInt_t          Vertex_fUniqueID[kMaxVertex];   //[Vertex_]
   UInt_t          Vertex_fBits[kMaxVertex];   //[Vertex_]
   Float_t         Vertex_T[kMaxVertex];   //[Vertex_]
   Float_t         Vertex_X[kMaxVertex];   //[Vertex_]
   Float_t         Vertex_Y[kMaxVertex];   //[Vertex_]
   Float_t         Vertex_Z[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_ErrorT[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_ErrorX[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_ErrorY[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_ErrorZ[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_Index[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_NDF[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_Sigma[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_SumPT2[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_GenSumPT2[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_GenDeltaZ[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_BTVSumPT2[kMaxVertex];   //[Vertex_]
   TRefArray       Vertex_Constituents[kMaxVertex];
   Int_t           Vertex_size;
   Int_t           GenJet_;
   UInt_t          GenJet_fUniqueID[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_fBits[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_PT[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Eta[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Phi[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_T[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Mass[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_DeltaEta[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_DeltaPhi[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_Flavor[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_FlavorAlgo[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_FlavorPhys[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_BTag[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_BTagAlgo[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_BTagPhys[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_TauTag[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_TauWeight[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_Charge[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_EhadOverEem[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NCharged[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NNeutrals[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Beta[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_BetaStar[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_MeanSqDeltaR[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_PTD[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_FracPt[kMaxGenJet][5];   //[GenJet_]
   Float_t         GenJet_Tau[kMaxGenJet][5];   //[GenJet_]
   TLorentzVector  GenJet_SoftDroppedJet[kMaxGenJet];
   TLorentzVector  GenJet_SoftDroppedSubJet1[kMaxGenJet];
   TLorentzVector  GenJet_SoftDroppedSubJet2[kMaxGenJet];
   TLorentzVector  GenJet_TrimmedP4[5][kMaxGenJet];
   TLorentzVector  GenJet_PrunedP4[5][kMaxGenJet];
   TLorentzVector  GenJet_SoftDroppedP4[5][kMaxGenJet];
   Int_t           GenJet_NSubJetsTrimmed[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NSubJetsPruned[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NSubJetsSoftDropped[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge23[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge34[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge45[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge56[kMaxGenJet];   //[GenJet_]
   TRefArray       GenJet_Constituents[kMaxGenJet];
   TRefArray       GenJet_Particles[kMaxGenJet];
   TLorentzVector  GenJet_Area[kMaxGenJet];
   Int_t           GenJet_size;
   Int_t           GenJetAK8_;
   UInt_t          GenJetAK8_fUniqueID[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_fBits[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_PT[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_Eta[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_Phi[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_T[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_Mass[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_DeltaEta[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_DeltaPhi[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_Flavor[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_FlavorAlgo[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_FlavorPhys[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_BTag[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_BTagAlgo[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_BTagPhys[kMaxGenJetAK8];   //[GenJetAK8_]
   UInt_t          GenJetAK8_TauTag[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_TauWeight[kMaxGenJetAK8];   //[GenJetAK8_]
   Int_t           GenJetAK8_Charge[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_EhadOverEem[kMaxGenJetAK8];   //[GenJetAK8_]
   Int_t           GenJetAK8_NCharged[kMaxGenJetAK8];   //[GenJetAK8_]
   Int_t           GenJetAK8_NNeutrals[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_Beta[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_BetaStar[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_MeanSqDeltaR[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_PTD[kMaxGenJetAK8];   //[GenJetAK8_]
   Float_t         GenJetAK8_FracPt[kMaxGenJetAK8][5];   //[GenJetAK8_]
   Float_t         GenJetAK8_Tau[kMaxGenJetAK8][5];   //[GenJetAK8_]
   TLorentzVector  GenJetAK8_SoftDroppedJet[kMaxGenJetAK8];
   TLorentzVector  GenJetAK8_SoftDroppedSubJet1[kMaxGenJetAK8];
   TLorentzVector  GenJetAK8_SoftDroppedSubJet2[kMaxGenJetAK8];
   TLorentzVector  GenJetAK8_TrimmedP4[5][kMaxGenJetAK8];
   TLorentzVector  GenJetAK8_PrunedP4[5][kMaxGenJetAK8];
   TLorentzVector  GenJetAK8_SoftDroppedP4[5][kMaxGenJetAK8];
   Int_t           GenJetAK8_NSubJetsTrimmed[kMaxGenJetAK8];   //[GenJetAK8_]
   Int_t           GenJetAK8_NSubJetsPruned[kMaxGenJetAK8];   //[GenJetAK8_]
   Int_t           GenJetAK8_NSubJetsSoftDropped[kMaxGenJetAK8];   //[GenJetAK8_]
   Double_t        GenJetAK8_ExclYmerge23[kMaxGenJetAK8];   //[GenJetAK8_]
   Double_t        GenJetAK8_ExclYmerge34[kMaxGenJetAK8];   //[GenJetAK8_]
   Double_t        GenJetAK8_ExclYmerge45[kMaxGenJetAK8];   //[GenJetAK8_]
   Double_t        GenJetAK8_ExclYmerge56[kMaxGenJetAK8];   //[GenJetAK8_]
   TRefArray       GenJetAK8_Constituents[kMaxGenJetAK8];
   TRefArray       GenJetAK8_Particles[kMaxGenJetAK8];
   TLorentzVector  GenJetAK8_Area[kMaxGenJetAK8];
   Int_t           GenJetAK8_size;
   Int_t           GenMissingET_;
   UInt_t          GenMissingET_fUniqueID[kMaxGenMissingET];   //[GenMissingET_]
   UInt_t          GenMissingET_fBits[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_MET[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_Eta[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_Phi[kMaxGenMissingET];   //[GenMissingET_]
   Int_t           GenMissingET_size;
   Int_t           PhotonLoose_;
   UInt_t          PhotonLoose_fUniqueID[kMaxPhotonLoose];   //[PhotonLoose_]
   UInt_t          PhotonLoose_fBits[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_PT[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_Eta[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_Phi[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_E[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_T[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_EhadOverEem[kMaxPhotonLoose];   //[PhotonLoose_]
   TRefArray       PhotonLoose_Particles[kMaxPhotonLoose];
   Float_t         PhotonLoose_IsolationVar[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_IsolationVarRhoCorr[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_SumPtCharged[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_SumPtNeutral[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_SumPtChargedPU[kMaxPhotonLoose];   //[PhotonLoose_]
   Float_t         PhotonLoose_SumPt[kMaxPhotonLoose];   //[PhotonLoose_]
   Int_t           PhotonLoose_Status[kMaxPhotonLoose];   //[PhotonLoose_]
   Int_t           PhotonLoose_size;
   Int_t           PhotonTight_;
   UInt_t          PhotonTight_fUniqueID[kMaxPhotonTight];   //[PhotonTight_]
   UInt_t          PhotonTight_fBits[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_PT[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_Eta[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_Phi[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_E[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_T[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_EhadOverEem[kMaxPhotonTight];   //[PhotonTight_]
   TRefArray       PhotonTight_Particles[kMaxPhotonTight];
   Float_t         PhotonTight_IsolationVar[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_IsolationVarRhoCorr[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_SumPtCharged[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_SumPtNeutral[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_SumPtChargedPU[kMaxPhotonTight];   //[PhotonTight_]
   Float_t         PhotonTight_SumPt[kMaxPhotonTight];   //[PhotonTight_]
   Int_t           PhotonTight_Status[kMaxPhotonTight];   //[PhotonTight_]
   Int_t           PhotonTight_size;
   Int_t           Electron_;
   UInt_t          Electron_fUniqueID[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fBits[kMaxElectron];   //[Electron_]
   Float_t         Electron_PT[kMaxElectron];   //[Electron_]
   Float_t         Electron_Eta[kMaxElectron];   //[Electron_]
   Float_t         Electron_Phi[kMaxElectron];   //[Electron_]
   Float_t         Electron_T[kMaxElectron];   //[Electron_]
   Int_t           Electron_Charge[kMaxElectron];   //[Electron_]
   Float_t         Electron_EhadOverEem[kMaxElectron];   //[Electron_]
   TRef            Electron_Particle[kMaxElectron];
   Float_t         Electron_IsolationVar[kMaxElectron];   //[Electron_]
   Float_t         Electron_IsolationVarRhoCorr[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtCharged[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtNeutral[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtChargedPU[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPt[kMaxElectron];   //[Electron_]
   Int_t           Electron_size;
   Int_t           MuonLoose_;
   UInt_t          MuonLoose_fUniqueID[kMaxMuonLoose];   //[MuonLoose_]
   UInt_t          MuonLoose_fBits[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_PT[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_Eta[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_Phi[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_T[kMaxMuonLoose];   //[MuonLoose_]
   Int_t           MuonLoose_Charge[kMaxMuonLoose];   //[MuonLoose_]
   TRef            MuonLoose_Particle[kMaxMuonLoose];
   Float_t         MuonLoose_IsolationVar[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_IsolationVarRhoCorr[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_SumPtCharged[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_SumPtNeutral[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_SumPtChargedPU[kMaxMuonLoose];   //[MuonLoose_]
   Float_t         MuonLoose_SumPt[kMaxMuonLoose];   //[MuonLoose_]
   Int_t           MuonLoose_size;
   Int_t           MuonTight_;
   UInt_t          MuonTight_fUniqueID[kMaxMuonTight];   //[MuonTight_]
   UInt_t          MuonTight_fBits[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_PT[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_Eta[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_Phi[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_T[kMaxMuonTight];   //[MuonTight_]
   Int_t           MuonTight_Charge[kMaxMuonTight];   //[MuonTight_]
   TRef            MuonTight_Particle[kMaxMuonTight];
   Float_t         MuonTight_IsolationVar[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_IsolationVarRhoCorr[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_SumPtCharged[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_SumPtNeutral[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_SumPtChargedPU[kMaxMuonTight];   //[MuonTight_]
   Float_t         MuonTight_SumPt[kMaxMuonTight];   //[MuonTight_]
   Int_t           MuonTight_size;
   Int_t           ElectronCHS_;
   UInt_t          ElectronCHS_fUniqueID[kMaxElectronCHS];   //[ElectronCHS_]
   UInt_t          ElectronCHS_fBits[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_PT[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_Eta[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_Phi[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_T[kMaxElectronCHS];   //[ElectronCHS_]
   Int_t           ElectronCHS_Charge[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_EhadOverEem[kMaxElectronCHS];   //[ElectronCHS_]
   TRef            ElectronCHS_Particle[kMaxElectronCHS];
   Float_t         ElectronCHS_IsolationVar[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_IsolationVarRhoCorr[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_SumPtCharged[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_SumPtNeutral[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_SumPtChargedPU[kMaxElectronCHS];   //[ElectronCHS_]
   Float_t         ElectronCHS_SumPt[kMaxElectronCHS];   //[ElectronCHS_]
   Int_t           ElectronCHS_size;
   Int_t           MuonLooseCHS_;
   UInt_t          MuonLooseCHS_fUniqueID[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   UInt_t          MuonLooseCHS_fBits[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_PT[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_Eta[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_Phi[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_T[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Int_t           MuonLooseCHS_Charge[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   TRef            MuonLooseCHS_Particle[kMaxMuonLooseCHS];
   Float_t         MuonLooseCHS_IsolationVar[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_IsolationVarRhoCorr[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_SumPtCharged[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_SumPtNeutral[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_SumPtChargedPU[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Float_t         MuonLooseCHS_SumPt[kMaxMuonLooseCHS];   //[MuonLooseCHS_]
   Int_t           MuonLooseCHS_size;
   Int_t           MuonTightCHS_;
   UInt_t          MuonTightCHS_fUniqueID[kMaxMuonTightCHS];   //[MuonTightCHS_]
   UInt_t          MuonTightCHS_fBits[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_PT[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_Eta[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_Phi[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_T[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Int_t           MuonTightCHS_Charge[kMaxMuonTightCHS];   //[MuonTightCHS_]
   TRef            MuonTightCHS_Particle[kMaxMuonTightCHS];
   Float_t         MuonTightCHS_IsolationVar[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_IsolationVarRhoCorr[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_SumPtCharged[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_SumPtNeutral[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_SumPtChargedPU[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Float_t         MuonTightCHS_SumPt[kMaxMuonTightCHS];   //[MuonTightCHS_]
   Int_t           MuonTightCHS_size;
   Int_t           Jet_;
   UInt_t          Jet_fUniqueID[kMaxJet];   //[Jet_]
   UInt_t          Jet_fBits[kMaxJet];   //[Jet_]
   Float_t         Jet_PT[kMaxJet];   //[Jet_]
   Float_t         Jet_Eta[kMaxJet];   //[Jet_]
   Float_t         Jet_Phi[kMaxJet];   //[Jet_]
   Float_t         Jet_T[kMaxJet];   //[Jet_]
   Float_t         Jet_Mass[kMaxJet];   //[Jet_]
   Float_t         Jet_DeltaEta[kMaxJet];   //[Jet_]
   Float_t         Jet_DeltaPhi[kMaxJet];   //[Jet_]
   UInt_t          Jet_Flavor[kMaxJet];   //[Jet_]
   UInt_t          Jet_FlavorAlgo[kMaxJet];   //[Jet_]
   UInt_t          Jet_FlavorPhys[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTag[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTagAlgo[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTagPhys[kMaxJet];   //[Jet_]
   UInt_t          Jet_TauTag[kMaxJet];   //[Jet_]
   Float_t         Jet_TauWeight[kMaxJet];   //[Jet_]
   Int_t           Jet_Charge[kMaxJet];   //[Jet_]
   Float_t         Jet_EhadOverEem[kMaxJet];   //[Jet_]
   Int_t           Jet_NCharged[kMaxJet];   //[Jet_]
   Int_t           Jet_NNeutrals[kMaxJet];   //[Jet_]
   Float_t         Jet_Beta[kMaxJet];   //[Jet_]
   Float_t         Jet_BetaStar[kMaxJet];   //[Jet_]
   Float_t         Jet_MeanSqDeltaR[kMaxJet];   //[Jet_]
   Float_t         Jet_PTD[kMaxJet];   //[Jet_]
   Float_t         Jet_FracPt[kMaxJet][5];   //[Jet_]
   Float_t         Jet_Tau[kMaxJet][5];   //[Jet_]
   TLorentzVector  Jet_SoftDroppedJet[kMaxJet];
   TLorentzVector  Jet_SoftDroppedSubJet1[kMaxJet];
   TLorentzVector  Jet_SoftDroppedSubJet2[kMaxJet];
   TLorentzVector  Jet_TrimmedP4[5][kMaxJet];
   TLorentzVector  Jet_PrunedP4[5][kMaxJet];
   TLorentzVector  Jet_SoftDroppedP4[5][kMaxJet];
   Int_t           Jet_NSubJetsTrimmed[kMaxJet];   //[Jet_]
   Int_t           Jet_NSubJetsPruned[kMaxJet];   //[Jet_]
   Int_t           Jet_NSubJetsSoftDropped[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge23[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge34[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge45[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge56[kMaxJet];   //[Jet_]
   TRefArray       Jet_Constituents[kMaxJet];
   TRefArray       Jet_Particles[kMaxJet];
   TLorentzVector  Jet_Area[kMaxJet];
   Int_t           Jet_size;
   Int_t           JetPUPPI_;
   UInt_t          JetPUPPI_fUniqueID[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_fBits[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_PT[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_Eta[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_Phi[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_T[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_Mass[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_DeltaEta[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_DeltaPhi[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_Flavor[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_FlavorAlgo[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_FlavorPhys[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_BTag[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_BTagAlgo[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_BTagPhys[kMaxJetPUPPI];   //[JetPUPPI_]
   UInt_t          JetPUPPI_TauTag[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_TauWeight[kMaxJetPUPPI];   //[JetPUPPI_]
   Int_t           JetPUPPI_Charge[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_EhadOverEem[kMaxJetPUPPI];   //[JetPUPPI_]
   Int_t           JetPUPPI_NCharged[kMaxJetPUPPI];   //[JetPUPPI_]
   Int_t           JetPUPPI_NNeutrals[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_Beta[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_BetaStar[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_MeanSqDeltaR[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_PTD[kMaxJetPUPPI];   //[JetPUPPI_]
   Float_t         JetPUPPI_FracPt[kMaxJetPUPPI][5];   //[JetPUPPI_]
   Float_t         JetPUPPI_Tau[kMaxJetPUPPI][5];   //[JetPUPPI_]
   TLorentzVector  JetPUPPI_SoftDroppedJet[kMaxJetPUPPI];
   TLorentzVector  JetPUPPI_SoftDroppedSubJet1[kMaxJetPUPPI];
   TLorentzVector  JetPUPPI_SoftDroppedSubJet2[kMaxJetPUPPI];
   TLorentzVector  JetPUPPI_TrimmedP4[5][kMaxJetPUPPI];
   TLorentzVector  JetPUPPI_PrunedP4[5][kMaxJetPUPPI];
   TLorentzVector  JetPUPPI_SoftDroppedP4[5][kMaxJetPUPPI];
   Int_t           JetPUPPI_NSubJetsTrimmed[kMaxJetPUPPI];   //[JetPUPPI_]
   Int_t           JetPUPPI_NSubJetsPruned[kMaxJetPUPPI];   //[JetPUPPI_]
   Int_t           JetPUPPI_NSubJetsSoftDropped[kMaxJetPUPPI];   //[JetPUPPI_]
   Double_t        JetPUPPI_ExclYmerge23[kMaxJetPUPPI];   //[JetPUPPI_]
   Double_t        JetPUPPI_ExclYmerge34[kMaxJetPUPPI];   //[JetPUPPI_]
   Double_t        JetPUPPI_ExclYmerge45[kMaxJetPUPPI];   //[JetPUPPI_]
   Double_t        JetPUPPI_ExclYmerge56[kMaxJetPUPPI];   //[JetPUPPI_]
   TRefArray       JetPUPPI_Constituents[kMaxJetPUPPI];
   TRefArray       JetPUPPI_Particles[kMaxJetPUPPI];
   TLorentzVector  JetPUPPI_Area[kMaxJetPUPPI];
   Int_t           JetPUPPI_size;
   Int_t           JetAK8_;
   UInt_t          JetAK8_fUniqueID[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_fBits[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_PT[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_Eta[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_Phi[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_T[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_Mass[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_DeltaEta[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_DeltaPhi[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_Flavor[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_FlavorAlgo[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_FlavorPhys[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_BTag[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_BTagAlgo[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_BTagPhys[kMaxJetAK8];   //[JetAK8_]
   UInt_t          JetAK8_TauTag[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_TauWeight[kMaxJetAK8];   //[JetAK8_]
   Int_t           JetAK8_Charge[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_EhadOverEem[kMaxJetAK8];   //[JetAK8_]
   Int_t           JetAK8_NCharged[kMaxJetAK8];   //[JetAK8_]
   Int_t           JetAK8_NNeutrals[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_Beta[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_BetaStar[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_MeanSqDeltaR[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_PTD[kMaxJetAK8];   //[JetAK8_]
   Float_t         JetAK8_FracPt[kMaxJetAK8][5];   //[JetAK8_]
   Float_t         JetAK8_Tau[kMaxJetAK8][5];   //[JetAK8_]
   TLorentzVector  JetAK8_SoftDroppedJet[kMaxJetAK8];
   TLorentzVector  JetAK8_SoftDroppedSubJet1[kMaxJetAK8];
   TLorentzVector  JetAK8_SoftDroppedSubJet2[kMaxJetAK8];
   TLorentzVector  JetAK8_TrimmedP4[5][kMaxJetAK8];
   TLorentzVector  JetAK8_PrunedP4[5][kMaxJetAK8];
   TLorentzVector  JetAK8_SoftDroppedP4[5][kMaxJetAK8];
   Int_t           JetAK8_NSubJetsTrimmed[kMaxJetAK8];   //[JetAK8_]
   Int_t           JetAK8_NSubJetsPruned[kMaxJetAK8];   //[JetAK8_]
   Int_t           JetAK8_NSubJetsSoftDropped[kMaxJetAK8];   //[JetAK8_]
   Double_t        JetAK8_ExclYmerge23[kMaxJetAK8];   //[JetAK8_]
   Double_t        JetAK8_ExclYmerge34[kMaxJetAK8];   //[JetAK8_]
   Double_t        JetAK8_ExclYmerge45[kMaxJetAK8];   //[JetAK8_]
   Double_t        JetAK8_ExclYmerge56[kMaxJetAK8];   //[JetAK8_]
   TRefArray       JetAK8_Constituents[kMaxJetAK8];
   TRefArray       JetAK8_Particles[kMaxJetAK8];
   TLorentzVector  JetAK8_Area[kMaxJetAK8];
   Int_t           JetAK8_size;
   Int_t           JetPUPPIAK8_;
   UInt_t          JetPUPPIAK8_fUniqueID[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_fBits[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_PT[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_Eta[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_Phi[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_T[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_Mass[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_DeltaEta[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_DeltaPhi[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_Flavor[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_FlavorAlgo[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_FlavorPhys[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_BTag[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_BTagAlgo[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_BTagPhys[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   UInt_t          JetPUPPIAK8_TauTag[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_TauWeight[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Int_t           JetPUPPIAK8_Charge[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_EhadOverEem[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Int_t           JetPUPPIAK8_NCharged[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Int_t           JetPUPPIAK8_NNeutrals[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_Beta[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_BetaStar[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_MeanSqDeltaR[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_PTD[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_FracPt[kMaxJetPUPPIAK8][5];   //[JetPUPPIAK8_]
   Float_t         JetPUPPIAK8_Tau[kMaxJetPUPPIAK8][5];   //[JetPUPPIAK8_]
   TLorentzVector  JetPUPPIAK8_SoftDroppedJet[kMaxJetPUPPIAK8];
   TLorentzVector  JetPUPPIAK8_SoftDroppedSubJet1[kMaxJetPUPPIAK8];
   TLorentzVector  JetPUPPIAK8_SoftDroppedSubJet2[kMaxJetPUPPIAK8];
   TLorentzVector  JetPUPPIAK8_TrimmedP4[5][kMaxJetPUPPIAK8];
   TLorentzVector  JetPUPPIAK8_PrunedP4[5][kMaxJetPUPPIAK8];
   TLorentzVector  JetPUPPIAK8_SoftDroppedP4[5][kMaxJetPUPPIAK8];
   Int_t           JetPUPPIAK8_NSubJetsTrimmed[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Int_t           JetPUPPIAK8_NSubJetsPruned[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Int_t           JetPUPPIAK8_NSubJetsSoftDropped[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Double_t        JetPUPPIAK8_ExclYmerge23[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Double_t        JetPUPPIAK8_ExclYmerge34[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Double_t        JetPUPPIAK8_ExclYmerge45[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   Double_t        JetPUPPIAK8_ExclYmerge56[kMaxJetPUPPIAK8];   //[JetPUPPIAK8_]
   TRefArray       JetPUPPIAK8_Constituents[kMaxJetPUPPIAK8];
   TRefArray       JetPUPPIAK8_Particles[kMaxJetPUPPIAK8];
   TLorentzVector  JetPUPPIAK8_Area[kMaxJetPUPPIAK8];
   Int_t           JetPUPPIAK8_size;
   Int_t           Rho_;
   UInt_t          Rho_fUniqueID[kMaxRho];   //[Rho_]
   UInt_t          Rho_fBits[kMaxRho];   //[Rho_]
   Float_t         Rho_Rho[kMaxRho];   //[Rho_]
   Float_t         Rho_Edges[kMaxRho][2];   //[Rho_]
   Int_t           Rho_size;
   Int_t           MissingET_;
   UInt_t          MissingET_fUniqueID[kMaxMissingET];   //[MissingET_]
   UInt_t          MissingET_fBits[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_MET[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_Eta[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_Phi[kMaxMissingET];   //[MissingET_]
   Int_t           MissingET_size;
   Int_t           PuppiMissingET_;
   UInt_t          PuppiMissingET_fUniqueID[kMaxPuppiMissingET];   //[PuppiMissingET_]
   UInt_t          PuppiMissingET_fBits[kMaxPuppiMissingET];   //[PuppiMissingET_]
   Float_t         PuppiMissingET_MET[kMaxPuppiMissingET];   //[PuppiMissingET_]
   Float_t         PuppiMissingET_Eta[kMaxPuppiMissingET];   //[PuppiMissingET_]
   Float_t         PuppiMissingET_Phi[kMaxPuppiMissingET];   //[PuppiMissingET_]
   Int_t           PuppiMissingET_size;
   Int_t           GenPileUpMissingET_;
   UInt_t          GenPileUpMissingET_fUniqueID[kMaxGenPileUpMissingET];   //[GenPileUpMissingET_]
   UInt_t          GenPileUpMissingET_fBits[kMaxGenPileUpMissingET];   //[GenPileUpMissingET_]
   Float_t         GenPileUpMissingET_MET[kMaxGenPileUpMissingET];   //[GenPileUpMissingET_]
   Float_t         GenPileUpMissingET_Eta[kMaxGenPileUpMissingET];   //[GenPileUpMissingET_]
   Float_t         GenPileUpMissingET_Phi[kMaxGenPileUpMissingET];   //[GenPileUpMissingET_]
   Int_t           GenPileUpMissingET_size;
   Int_t           ScalarHT_;
   UInt_t          ScalarHT_fUniqueID[kMaxScalarHT];   //[ScalarHT_]
   UInt_t          ScalarHT_fBits[kMaxScalarHT];   //[ScalarHT_]
   Float_t         ScalarHT_HT[kMaxScalarHT];   //[ScalarHT_]
   Int_t           ScalarHT_size;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_ReadTime;   //!
   TBranch        *b_Event_ProcTime;   //!
   TBranch        *b_Event_ProcessID;   //!
   TBranch        *b_Event_MPI;   //!
   TBranch        *b_Event_Weight;   //!
   TBranch        *b_Event_Scale;   //!
   TBranch        *b_Event_AlphaQED;   //!
   TBranch        *b_Event_AlphaQCD;   //!
   TBranch        *b_Event_ID1;   //!
   TBranch        *b_Event_ID2;   //!
   TBranch        *b_Event_X1;   //!
   TBranch        *b_Event_X2;   //!
   TBranch        *b_Event_ScalePDF;   //!
   TBranch        *b_Event_PDF1;   //!
   TBranch        *b_Event_PDF2;   //!
   TBranch        *b_Event_size;   //!
   TBranch        *b_Weight_;   //!
   TBranch        *b_Weight_fUniqueID;   //!
   TBranch        *b_Weight_fBits;   //!
   TBranch        *b_Weight_Weight;   //!
   TBranch        *b_Weight_size;   //!
   TBranch        *b_Particle_;   //!
   TBranch        *b_Particle_fUniqueID;   //!
   TBranch        *b_Particle_fBits;   //!
   TBranch        *b_Particle_PID;   //!
   TBranch        *b_Particle_Status;   //!
   TBranch        *b_Particle_IsPU;   //!
   TBranch        *b_Particle_M1;   //!
   TBranch        *b_Particle_M2;   //!
   TBranch        *b_Particle_D1;   //!
   TBranch        *b_Particle_D2;   //!
   TBranch        *b_Particle_Charge;   //!
   TBranch        *b_Particle_Mass;   //!
   TBranch        *b_Particle_E;   //!
   TBranch        *b_Particle_Px;   //!
   TBranch        *b_Particle_Py;   //!
   TBranch        *b_Particle_Pz;   //!
   TBranch        *b_Particle_P;   //!
   TBranch        *b_Particle_PT;   //!
   TBranch        *b_Particle_Eta;   //!
   TBranch        *b_Particle_Phi;   //!
   TBranch        *b_Particle_Rapidity;   //!
   TBranch        *b_Particle_CtgTheta;   //!
   TBranch        *b_Particle_D0;   //!
   TBranch        *b_Particle_DZ;   //!
   TBranch        *b_Particle_T;   //!
   TBranch        *b_Particle_X;   //!
   TBranch        *b_Particle_Y;   //!
   TBranch        *b_Particle_Z;   //!
   TBranch        *b_Particle_size;   //!
   TBranch        *b_Vertex_;   //!
   TBranch        *b_Vertex_fUniqueID;   //!
   TBranch        *b_Vertex_fBits;   //!
   TBranch        *b_Vertex_T;   //!
   TBranch        *b_Vertex_X;   //!
   TBranch        *b_Vertex_Y;   //!
   TBranch        *b_Vertex_Z;   //!
   TBranch        *b_Vertex_ErrorT;   //!
   TBranch        *b_Vertex_ErrorX;   //!
   TBranch        *b_Vertex_ErrorY;   //!
   TBranch        *b_Vertex_ErrorZ;   //!
   TBranch        *b_Vertex_Index;   //!
   TBranch        *b_Vertex_NDF;   //!
   TBranch        *b_Vertex_Sigma;   //!
   TBranch        *b_Vertex_SumPT2;   //!
   TBranch        *b_Vertex_GenSumPT2;   //!
   TBranch        *b_Vertex_GenDeltaZ;   //!
   TBranch        *b_Vertex_BTVSumPT2;   //!
   TBranch        *b_Vertex_Constituents;   //!
   TBranch        *b_Vertex_size;   //!
   TBranch        *b_GenJet_;   //!
   TBranch        *b_GenJet_fUniqueID;   //!
   TBranch        *b_GenJet_fBits;   //!
   TBranch        *b_GenJet_PT;   //!
   TBranch        *b_GenJet_Eta;   //!
   TBranch        *b_GenJet_Phi;   //!
   TBranch        *b_GenJet_T;   //!
   TBranch        *b_GenJet_Mass;   //!
   TBranch        *b_GenJet_DeltaEta;   //!
   TBranch        *b_GenJet_DeltaPhi;   //!
   TBranch        *b_GenJet_Flavor;   //!
   TBranch        *b_GenJet_FlavorAlgo;   //!
   TBranch        *b_GenJet_FlavorPhys;   //!
   TBranch        *b_GenJet_BTag;   //!
   TBranch        *b_GenJet_BTagAlgo;   //!
   TBranch        *b_GenJet_BTagPhys;   //!
   TBranch        *b_GenJet_TauTag;   //!
   TBranch        *b_GenJet_TauWeight;   //!
   TBranch        *b_GenJet_Charge;   //!
   TBranch        *b_GenJet_EhadOverEem;   //!
   TBranch        *b_GenJet_NCharged;   //!
   TBranch        *b_GenJet_NNeutrals;   //!
   TBranch        *b_GenJet_Beta;   //!
   TBranch        *b_GenJet_BetaStar;   //!
   TBranch        *b_GenJet_MeanSqDeltaR;   //!
   TBranch        *b_GenJet_PTD;   //!
   TBranch        *b_GenJet_FracPt;   //!
   TBranch        *b_GenJet_Tau;   //!
   TBranch        *b_GenJet_SoftDroppedJet;   //!
   TBranch        *b_GenJet_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJet_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJet_TrimmedP4;   //!
   TBranch        *b_GenJet_PrunedP4;   //!
   TBranch        *b_GenJet_SoftDroppedP4;   //!
   TBranch        *b_GenJet_NSubJetsTrimmed;   //!
   TBranch        *b_GenJet_NSubJetsPruned;   //!
   TBranch        *b_GenJet_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJet_ExclYmerge23;   //!
   TBranch        *b_GenJet_ExclYmerge34;   //!
   TBranch        *b_GenJet_ExclYmerge45;   //!
   TBranch        *b_GenJet_ExclYmerge56;   //!
   TBranch        *b_GenJet_Constituents;   //!
   TBranch        *b_GenJet_Particles;   //!
   TBranch        *b_GenJet_Area;   //!
   TBranch        *b_GenJet_size;   //!
   TBranch        *b_GenJetAK8_;   //!
   TBranch        *b_GenJetAK8_fUniqueID;   //!
   TBranch        *b_GenJetAK8_fBits;   //!
   TBranch        *b_GenJetAK8_PT;   //!
   TBranch        *b_GenJetAK8_Eta;   //!
   TBranch        *b_GenJetAK8_Phi;   //!
   TBranch        *b_GenJetAK8_T;   //!
   TBranch        *b_GenJetAK8_Mass;   //!
   TBranch        *b_GenJetAK8_DeltaEta;   //!
   TBranch        *b_GenJetAK8_DeltaPhi;   //!
   TBranch        *b_GenJetAK8_Flavor;   //!
   TBranch        *b_GenJetAK8_FlavorAlgo;   //!
   TBranch        *b_GenJetAK8_FlavorPhys;   //!
   TBranch        *b_GenJetAK8_BTag;   //!
   TBranch        *b_GenJetAK8_BTagAlgo;   //!
   TBranch        *b_GenJetAK8_BTagPhys;   //!
   TBranch        *b_GenJetAK8_TauTag;   //!
   TBranch        *b_GenJetAK8_TauWeight;   //!
   TBranch        *b_GenJetAK8_Charge;   //!
   TBranch        *b_GenJetAK8_EhadOverEem;   //!
   TBranch        *b_GenJetAK8_NCharged;   //!
   TBranch        *b_GenJetAK8_NNeutrals;   //!
   TBranch        *b_GenJetAK8_Beta;   //!
   TBranch        *b_GenJetAK8_BetaStar;   //!
   TBranch        *b_GenJetAK8_MeanSqDeltaR;   //!
   TBranch        *b_GenJetAK8_PTD;   //!
   TBranch        *b_GenJetAK8_FracPt;   //!
   TBranch        *b_GenJetAK8_Tau;   //!
   TBranch        *b_GenJetAK8_SoftDroppedJet;   //!
   TBranch        *b_GenJetAK8_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJetAK8_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJetAK8_TrimmedP4;   //!
   TBranch        *b_GenJetAK8_PrunedP4;   //!
   TBranch        *b_GenJetAK8_SoftDroppedP4;   //!
   TBranch        *b_GenJetAK8_NSubJetsTrimmed;   //!
   TBranch        *b_GenJetAK8_NSubJetsPruned;   //!
   TBranch        *b_GenJetAK8_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJetAK8_ExclYmerge23;   //!
   TBranch        *b_GenJetAK8_ExclYmerge34;   //!
   TBranch        *b_GenJetAK8_ExclYmerge45;   //!
   TBranch        *b_GenJetAK8_ExclYmerge56;   //!
   TBranch        *b_GenJetAK8_Constituents;   //!
   TBranch        *b_GenJetAK8_Particles;   //!
   TBranch        *b_GenJetAK8_Area;   //!
   TBranch        *b_GenJetAK8_size;   //!
   TBranch        *b_GenMissingET_;   //!
   TBranch        *b_GenMissingET_fUniqueID;   //!
   TBranch        *b_GenMissingET_fBits;   //!
   TBranch        *b_GenMissingET_MET;   //!
   TBranch        *b_GenMissingET_Eta;   //!
   TBranch        *b_GenMissingET_Phi;   //!
   TBranch        *b_GenMissingET_size;   //!
   TBranch        *b_PhotonLoose_;   //!
   TBranch        *b_PhotonLoose_fUniqueID;   //!
   TBranch        *b_PhotonLoose_fBits;   //!
   TBranch        *b_PhotonLoose_PT;   //!
   TBranch        *b_PhotonLoose_Eta;   //!
   TBranch        *b_PhotonLoose_Phi;   //!
   TBranch        *b_PhotonLoose_E;   //!
   TBranch        *b_PhotonLoose_T;   //!
   TBranch        *b_PhotonLoose_EhadOverEem;   //!
   TBranch        *b_PhotonLoose_Particles;   //!
   TBranch        *b_PhotonLoose_IsolationVar;   //!
   TBranch        *b_PhotonLoose_IsolationVarRhoCorr;   //!
   TBranch        *b_PhotonLoose_SumPtCharged;   //!
   TBranch        *b_PhotonLoose_SumPtNeutral;   //!
   TBranch        *b_PhotonLoose_SumPtChargedPU;   //!
   TBranch        *b_PhotonLoose_SumPt;   //!
   TBranch        *b_PhotonLoose_Status;   //!
   TBranch        *b_PhotonLoose_size;   //!
   TBranch        *b_PhotonTight_;   //!
   TBranch        *b_PhotonTight_fUniqueID;   //!
   TBranch        *b_PhotonTight_fBits;   //!
   TBranch        *b_PhotonTight_PT;   //!
   TBranch        *b_PhotonTight_Eta;   //!
   TBranch        *b_PhotonTight_Phi;   //!
   TBranch        *b_PhotonTight_E;   //!
   TBranch        *b_PhotonTight_T;   //!
   TBranch        *b_PhotonTight_EhadOverEem;   //!
   TBranch        *b_PhotonTight_Particles;   //!
   TBranch        *b_PhotonTight_IsolationVar;   //!
   TBranch        *b_PhotonTight_IsolationVarRhoCorr;   //!
   TBranch        *b_PhotonTight_SumPtCharged;   //!
   TBranch        *b_PhotonTight_SumPtNeutral;   //!
   TBranch        *b_PhotonTight_SumPtChargedPU;   //!
   TBranch        *b_PhotonTight_SumPt;   //!
   TBranch        *b_PhotonTight_Status;   //!
   TBranch        *b_PhotonTight_size;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_fUniqueID;   //!
   TBranch        *b_Electron_fBits;   //!
   TBranch        *b_Electron_PT;   //!
   TBranch        *b_Electron_Eta;   //!
   TBranch        *b_Electron_Phi;   //!
   TBranch        *b_Electron_T;   //!
   TBranch        *b_Electron_Charge;   //!
   TBranch        *b_Electron_EhadOverEem;   //!
   TBranch        *b_Electron_Particle;   //!
   TBranch        *b_Electron_IsolationVar;   //!
   TBranch        *b_Electron_IsolationVarRhoCorr;   //!
   TBranch        *b_Electron_SumPtCharged;   //!
   TBranch        *b_Electron_SumPtNeutral;   //!
   TBranch        *b_Electron_SumPtChargedPU;   //!
   TBranch        *b_Electron_SumPt;   //!
   TBranch        *b_Electron_size;   //!
   TBranch        *b_MuonLoose_;   //!
   TBranch        *b_MuonLoose_fUniqueID;   //!
   TBranch        *b_MuonLoose_fBits;   //!
   TBranch        *b_MuonLoose_PT;   //!
   TBranch        *b_MuonLoose_Eta;   //!
   TBranch        *b_MuonLoose_Phi;   //!
   TBranch        *b_MuonLoose_T;   //!
   TBranch        *b_MuonLoose_Charge;   //!
   TBranch        *b_MuonLoose_Particle;   //!
   TBranch        *b_MuonLoose_IsolationVar;   //!
   TBranch        *b_MuonLoose_IsolationVarRhoCorr;   //!
   TBranch        *b_MuonLoose_SumPtCharged;   //!
   TBranch        *b_MuonLoose_SumPtNeutral;   //!
   TBranch        *b_MuonLoose_SumPtChargedPU;   //!
   TBranch        *b_MuonLoose_SumPt;   //!
   TBranch        *b_MuonLoose_size;   //!
   TBranch        *b_MuonTight_;   //!
   TBranch        *b_MuonTight_fUniqueID;   //!
   TBranch        *b_MuonTight_fBits;   //!
   TBranch        *b_MuonTight_PT;   //!
   TBranch        *b_MuonTight_Eta;   //!
   TBranch        *b_MuonTight_Phi;   //!
   TBranch        *b_MuonTight_T;   //!
   TBranch        *b_MuonTight_Charge;   //!
   TBranch        *b_MuonTight_Particle;   //!
   TBranch        *b_MuonTight_IsolationVar;   //!
   TBranch        *b_MuonTight_IsolationVarRhoCorr;   //!
   TBranch        *b_MuonTight_SumPtCharged;   //!
   TBranch        *b_MuonTight_SumPtNeutral;   //!
   TBranch        *b_MuonTight_SumPtChargedPU;   //!
   TBranch        *b_MuonTight_SumPt;   //!
   TBranch        *b_MuonTight_size;   //!
   TBranch        *b_ElectronCHS_;   //!
   TBranch        *b_ElectronCHS_fUniqueID;   //!
   TBranch        *b_ElectronCHS_fBits;   //!
   TBranch        *b_ElectronCHS_PT;   //!
   TBranch        *b_ElectronCHS_Eta;   //!
   TBranch        *b_ElectronCHS_Phi;   //!
   TBranch        *b_ElectronCHS_T;   //!
   TBranch        *b_ElectronCHS_Charge;   //!
   TBranch        *b_ElectronCHS_EhadOverEem;   //!
   TBranch        *b_ElectronCHS_Particle;   //!
   TBranch        *b_ElectronCHS_IsolationVar;   //!
   TBranch        *b_ElectronCHS_IsolationVarRhoCorr;   //!
   TBranch        *b_ElectronCHS_SumPtCharged;   //!
   TBranch        *b_ElectronCHS_SumPtNeutral;   //!
   TBranch        *b_ElectronCHS_SumPtChargedPU;   //!
   TBranch        *b_ElectronCHS_SumPt;   //!
   TBranch        *b_ElectronCHS_size;   //!
   TBranch        *b_MuonLooseCHS_;   //!
   TBranch        *b_MuonLooseCHS_fUniqueID;   //!
   TBranch        *b_MuonLooseCHS_fBits;   //!
   TBranch        *b_MuonLooseCHS_PT;   //!
   TBranch        *b_MuonLooseCHS_Eta;   //!
   TBranch        *b_MuonLooseCHS_Phi;   //!
   TBranch        *b_MuonLooseCHS_T;   //!
   TBranch        *b_MuonLooseCHS_Charge;   //!
   TBranch        *b_MuonLooseCHS_Particle;   //!
   TBranch        *b_MuonLooseCHS_IsolationVar;   //!
   TBranch        *b_MuonLooseCHS_IsolationVarRhoCorr;   //!
   TBranch        *b_MuonLooseCHS_SumPtCharged;   //!
   TBranch        *b_MuonLooseCHS_SumPtNeutral;   //!
   TBranch        *b_MuonLooseCHS_SumPtChargedPU;   //!
   TBranch        *b_MuonLooseCHS_SumPt;   //!
   TBranch        *b_MuonLooseCHS_size;   //!
   TBranch        *b_MuonTightCHS_;   //!
   TBranch        *b_MuonTightCHS_fUniqueID;   //!
   TBranch        *b_MuonTightCHS_fBits;   //!
   TBranch        *b_MuonTightCHS_PT;   //!
   TBranch        *b_MuonTightCHS_Eta;   //!
   TBranch        *b_MuonTightCHS_Phi;   //!
   TBranch        *b_MuonTightCHS_T;   //!
   TBranch        *b_MuonTightCHS_Charge;   //!
   TBranch        *b_MuonTightCHS_Particle;   //!
   TBranch        *b_MuonTightCHS_IsolationVar;   //!
   TBranch        *b_MuonTightCHS_IsolationVarRhoCorr;   //!
   TBranch        *b_MuonTightCHS_SumPtCharged;   //!
   TBranch        *b_MuonTightCHS_SumPtNeutral;   //!
   TBranch        *b_MuonTightCHS_SumPtChargedPU;   //!
   TBranch        *b_MuonTightCHS_SumPt;   //!
   TBranch        *b_MuonTightCHS_size;   //!
   TBranch        *b_Jet_;   //!
   TBranch        *b_Jet_fUniqueID;   //!
   TBranch        *b_Jet_fBits;   //!
   TBranch        *b_Jet_PT;   //!
   TBranch        *b_Jet_Eta;   //!
   TBranch        *b_Jet_Phi;   //!
   TBranch        *b_Jet_T;   //!
   TBranch        *b_Jet_Mass;   //!
   TBranch        *b_Jet_DeltaEta;   //!
   TBranch        *b_Jet_DeltaPhi;   //!
   TBranch        *b_Jet_Flavor;   //!
   TBranch        *b_Jet_FlavorAlgo;   //!
   TBranch        *b_Jet_FlavorPhys;   //!
   TBranch        *b_Jet_BTag;   //!
   TBranch        *b_Jet_BTagAlgo;   //!
   TBranch        *b_Jet_BTagPhys;   //!
   TBranch        *b_Jet_TauTag;   //!
   TBranch        *b_Jet_TauWeight;   //!
   TBranch        *b_Jet_Charge;   //!
   TBranch        *b_Jet_EhadOverEem;   //!
   TBranch        *b_Jet_NCharged;   //!
   TBranch        *b_Jet_NNeutrals;   //!
   TBranch        *b_Jet_Beta;   //!
   TBranch        *b_Jet_BetaStar;   //!
   TBranch        *b_Jet_MeanSqDeltaR;   //!
   TBranch        *b_Jet_PTD;   //!
   TBranch        *b_Jet_FracPt;   //!
   TBranch        *b_Jet_Tau;   //!
   TBranch        *b_Jet_SoftDroppedJet;   //!
   TBranch        *b_Jet_SoftDroppedSubJet1;   //!
   TBranch        *b_Jet_SoftDroppedSubJet2;   //!
   TBranch        *b_Jet_TrimmedP4;   //!
   TBranch        *b_Jet_PrunedP4;   //!
   TBranch        *b_Jet_SoftDroppedP4;   //!
   TBranch        *b_Jet_NSubJetsTrimmed;   //!
   TBranch        *b_Jet_NSubJetsPruned;   //!
   TBranch        *b_Jet_NSubJetsSoftDropped;   //!
   TBranch        *b_Jet_ExclYmerge23;   //!
   TBranch        *b_Jet_ExclYmerge34;   //!
   TBranch        *b_Jet_ExclYmerge45;   //!
   TBranch        *b_Jet_ExclYmerge56;   //!
   TBranch        *b_Jet_Constituents;   //!
   TBranch        *b_Jet_Particles;   //!
   TBranch        *b_Jet_Area;   //!
   TBranch        *b_Jet_size;   //!
   TBranch        *b_JetPUPPI_;   //!
   TBranch        *b_JetPUPPI_fUniqueID;   //!
   TBranch        *b_JetPUPPI_fBits;   //!
   TBranch        *b_JetPUPPI_PT;   //!
   TBranch        *b_JetPUPPI_Eta;   //!
   TBranch        *b_JetPUPPI_Phi;   //!
   TBranch        *b_JetPUPPI_T;   //!
   TBranch        *b_JetPUPPI_Mass;   //!
   TBranch        *b_JetPUPPI_DeltaEta;   //!
   TBranch        *b_JetPUPPI_DeltaPhi;   //!
   TBranch        *b_JetPUPPI_Flavor;   //!
   TBranch        *b_JetPUPPI_FlavorAlgo;   //!
   TBranch        *b_JetPUPPI_FlavorPhys;   //!
   TBranch        *b_JetPUPPI_BTag;   //!
   TBranch        *b_JetPUPPI_BTagAlgo;   //!
   TBranch        *b_JetPUPPI_BTagPhys;   //!
   TBranch        *b_JetPUPPI_TauTag;   //!
   TBranch        *b_JetPUPPI_TauWeight;   //!
   TBranch        *b_JetPUPPI_Charge;   //!
   TBranch        *b_JetPUPPI_EhadOverEem;   //!
   TBranch        *b_JetPUPPI_NCharged;   //!
   TBranch        *b_JetPUPPI_NNeutrals;   //!
   TBranch        *b_JetPUPPI_Beta;   //!
   TBranch        *b_JetPUPPI_BetaStar;   //!
   TBranch        *b_JetPUPPI_MeanSqDeltaR;   //!
   TBranch        *b_JetPUPPI_PTD;   //!
   TBranch        *b_JetPUPPI_FracPt;   //!
   TBranch        *b_JetPUPPI_Tau;   //!
   TBranch        *b_JetPUPPI_SoftDroppedJet;   //!
   TBranch        *b_JetPUPPI_SoftDroppedSubJet1;   //!
   TBranch        *b_JetPUPPI_SoftDroppedSubJet2;   //!
   TBranch        *b_JetPUPPI_TrimmedP4;   //!
   TBranch        *b_JetPUPPI_PrunedP4;   //!
   TBranch        *b_JetPUPPI_SoftDroppedP4;   //!
   TBranch        *b_JetPUPPI_NSubJetsTrimmed;   //!
   TBranch        *b_JetPUPPI_NSubJetsPruned;   //!
   TBranch        *b_JetPUPPI_NSubJetsSoftDropped;   //!
   TBranch        *b_JetPUPPI_ExclYmerge23;   //!
   TBranch        *b_JetPUPPI_ExclYmerge34;   //!
   TBranch        *b_JetPUPPI_ExclYmerge45;   //!
   TBranch        *b_JetPUPPI_ExclYmerge56;   //!
   TBranch        *b_JetPUPPI_Constituents;   //!
   TBranch        *b_JetPUPPI_Particles;   //!
   TBranch        *b_JetPUPPI_Area;   //!
   TBranch        *b_JetPUPPI_size;   //!
   TBranch        *b_JetAK8_;   //!
   TBranch        *b_JetAK8_fUniqueID;   //!
   TBranch        *b_JetAK8_fBits;   //!
   TBranch        *b_JetAK8_PT;   //!
   TBranch        *b_JetAK8_Eta;   //!
   TBranch        *b_JetAK8_Phi;   //!
   TBranch        *b_JetAK8_T;   //!
   TBranch        *b_JetAK8_Mass;   //!
   TBranch        *b_JetAK8_DeltaEta;   //!
   TBranch        *b_JetAK8_DeltaPhi;   //!
   TBranch        *b_JetAK8_Flavor;   //!
   TBranch        *b_JetAK8_FlavorAlgo;   //!
   TBranch        *b_JetAK8_FlavorPhys;   //!
   TBranch        *b_JetAK8_BTag;   //!
   TBranch        *b_JetAK8_BTagAlgo;   //!
   TBranch        *b_JetAK8_BTagPhys;   //!
   TBranch        *b_JetAK8_TauTag;   //!
   TBranch        *b_JetAK8_TauWeight;   //!
   TBranch        *b_JetAK8_Charge;   //!
   TBranch        *b_JetAK8_EhadOverEem;   //!
   TBranch        *b_JetAK8_NCharged;   //!
   TBranch        *b_JetAK8_NNeutrals;   //!
   TBranch        *b_JetAK8_Beta;   //!
   TBranch        *b_JetAK8_BetaStar;   //!
   TBranch        *b_JetAK8_MeanSqDeltaR;   //!
   TBranch        *b_JetAK8_PTD;   //!
   TBranch        *b_JetAK8_FracPt;   //!
   TBranch        *b_JetAK8_Tau;   //!
   TBranch        *b_JetAK8_SoftDroppedJet;   //!
   TBranch        *b_JetAK8_SoftDroppedSubJet1;   //!
   TBranch        *b_JetAK8_SoftDroppedSubJet2;   //!
   TBranch        *b_JetAK8_TrimmedP4;   //!
   TBranch        *b_JetAK8_PrunedP4;   //!
   TBranch        *b_JetAK8_SoftDroppedP4;   //!
   TBranch        *b_JetAK8_NSubJetsTrimmed;   //!
   TBranch        *b_JetAK8_NSubJetsPruned;   //!
   TBranch        *b_JetAK8_NSubJetsSoftDropped;   //!
   TBranch        *b_JetAK8_ExclYmerge23;   //!
   TBranch        *b_JetAK8_ExclYmerge34;   //!
   TBranch        *b_JetAK8_ExclYmerge45;   //!
   TBranch        *b_JetAK8_ExclYmerge56;   //!
   TBranch        *b_JetAK8_Constituents;   //!
   TBranch        *b_JetAK8_Particles;   //!
   TBranch        *b_JetAK8_Area;   //!
   TBranch        *b_JetAK8_size;   //!
   TBranch        *b_JetPUPPIAK8_;   //!
   TBranch        *b_JetPUPPIAK8_fUniqueID;   //!
   TBranch        *b_JetPUPPIAK8_fBits;   //!
   TBranch        *b_JetPUPPIAK8_PT;   //!
   TBranch        *b_JetPUPPIAK8_Eta;   //!
   TBranch        *b_JetPUPPIAK8_Phi;   //!
   TBranch        *b_JetPUPPIAK8_T;   //!
   TBranch        *b_JetPUPPIAK8_Mass;   //!
   TBranch        *b_JetPUPPIAK8_DeltaEta;   //!
   TBranch        *b_JetPUPPIAK8_DeltaPhi;   //!
   TBranch        *b_JetPUPPIAK8_Flavor;   //!
   TBranch        *b_JetPUPPIAK8_FlavorAlgo;   //!
   TBranch        *b_JetPUPPIAK8_FlavorPhys;   //!
   TBranch        *b_JetPUPPIAK8_BTag;   //!
   TBranch        *b_JetPUPPIAK8_BTagAlgo;   //!
   TBranch        *b_JetPUPPIAK8_BTagPhys;   //!
   TBranch        *b_JetPUPPIAK8_TauTag;   //!
   TBranch        *b_JetPUPPIAK8_TauWeight;   //!
   TBranch        *b_JetPUPPIAK8_Charge;   //!
   TBranch        *b_JetPUPPIAK8_EhadOverEem;   //!
   TBranch        *b_JetPUPPIAK8_NCharged;   //!
   TBranch        *b_JetPUPPIAK8_NNeutrals;   //!
   TBranch        *b_JetPUPPIAK8_Beta;   //!
   TBranch        *b_JetPUPPIAK8_BetaStar;   //!
   TBranch        *b_JetPUPPIAK8_MeanSqDeltaR;   //!
   TBranch        *b_JetPUPPIAK8_PTD;   //!
   TBranch        *b_JetPUPPIAK8_FracPt;   //!
   TBranch        *b_JetPUPPIAK8_Tau;   //!
   TBranch        *b_JetPUPPIAK8_SoftDroppedJet;   //!
   TBranch        *b_JetPUPPIAK8_SoftDroppedSubJet1;   //!
   TBranch        *b_JetPUPPIAK8_SoftDroppedSubJet2;   //!
   TBranch        *b_JetPUPPIAK8_TrimmedP4;   //!
   TBranch        *b_JetPUPPIAK8_PrunedP4;   //!
   TBranch        *b_JetPUPPIAK8_SoftDroppedP4;   //!
   TBranch        *b_JetPUPPIAK8_NSubJetsTrimmed;   //!
   TBranch        *b_JetPUPPIAK8_NSubJetsPruned;   //!
   TBranch        *b_JetPUPPIAK8_NSubJetsSoftDropped;   //!
   TBranch        *b_JetPUPPIAK8_ExclYmerge23;   //!
   TBranch        *b_JetPUPPIAK8_ExclYmerge34;   //!
   TBranch        *b_JetPUPPIAK8_ExclYmerge45;   //!
   TBranch        *b_JetPUPPIAK8_ExclYmerge56;   //!
   TBranch        *b_JetPUPPIAK8_Constituents;   //!
   TBranch        *b_JetPUPPIAK8_Particles;   //!
   TBranch        *b_JetPUPPIAK8_Area;   //!
   TBranch        *b_JetPUPPIAK8_size;   //!
   TBranch        *b_Rho_;   //!
   TBranch        *b_Rho_fUniqueID;   //!
   TBranch        *b_Rho_fBits;   //!
   TBranch        *b_Rho_Rho;   //!
   TBranch        *b_Rho_Edges;   //!
   TBranch        *b_Rho_size;   //!
   TBranch        *b_MissingET_;   //!
   TBranch        *b_MissingET_fUniqueID;   //!
   TBranch        *b_MissingET_fBits;   //!
   TBranch        *b_MissingET_MET;   //!
   TBranch        *b_MissingET_Eta;   //!
   TBranch        *b_MissingET_Phi;   //!
   TBranch        *b_MissingET_size;   //!
   TBranch        *b_PuppiMissingET_;   //!
   TBranch        *b_PuppiMissingET_fUniqueID;   //!
   TBranch        *b_PuppiMissingET_fBits;   //!
   TBranch        *b_PuppiMissingET_MET;   //!
   TBranch        *b_PuppiMissingET_Eta;   //!
   TBranch        *b_PuppiMissingET_Phi;   //!
   TBranch        *b_PuppiMissingET_size;   //!
   TBranch        *b_GenPileUpMissingET_;   //!
   TBranch        *b_GenPileUpMissingET_fUniqueID;   //!
   TBranch        *b_GenPileUpMissingET_fBits;   //!
   TBranch        *b_GenPileUpMissingET_MET;   //!
   TBranch        *b_GenPileUpMissingET_Eta;   //!
   TBranch        *b_GenPileUpMissingET_Phi;   //!
   TBranch        *b_GenPileUpMissingET_size;   //!
   TBranch        *b_ScalarHT_;   //!
   TBranch        *b_ScalarHT_fUniqueID;   //!
   TBranch        *b_ScalarHT_fBits;   //!
   TBranch        *b_ScalarHT_HT;   //!
   TBranch        *b_ScalarHT_size;   //!

   delphesReader(TTree *tree=0);
   virtual ~delphesReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif