/*
 * main.cc
 *
 *  Created on: 8 Apr 2016
 *      Author: giles giles.strong@outlook.com
 */
//Local
#include "main.hh"

bool debug = false;

bool getOSTauLeptonPair(delphesReader* reader, std::vector<int>* taus,	std::vector<int> *leptons,
		int* tau, int* lepton, bool electron) {
	/*Checks whether an OS tau-lepton pair exists, and if so returns true and points tau and lepton
 	to the selected particles*/
	std::vector<std::pair<int, int> > pairs; //Initialise array for tau lepton OS pairs
	for (int t : *taus) { //Loop through taus
		for (int l : *leptons) { //Loop through leptons
			if (electron) {
				if (reader->Jet_Charge[t] != reader->Electron_Charge[l]) {
					pairs.push_back(std::make_pair(t, l));
				}
			} else {
				if (reader->Jet_Charge[t] != reader->Muon_Charge[l]) {
					pairs.push_back(std::make_pair(t, l));
				}
			}
		}
		if (pairs.size() == 1) { //Only one OS pair found
			*tau = pairs[0].first;
			*lepton = pairs[0].second;
			return true;
		} else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
			double pT, highestPT = 0;
			std::pair<int, int> best;
			for (std::pair<int, int> p : pairs) { //Loop through pairs
				if (electron) {
					pT = std::abs(reader->Jet_PT[p.first])+std::abs(reader->Electron_PT[p.second]);
				} else {
					pT = std::abs(reader->Jet_PT[p.first])+std::abs(reader->Muon_PT[p.second]);
				}
				if (pT > highestPT) {
					best = p;
					highestPT = pT;
				}
			}
			*tau = best.first;
			*lepton = best.second;
			return true;
		} else { //No OS pairs found
			return false;
		}
	}
	return false;
}

bool getOSLeptonLeptonPair(delphesReader* reader, std::vector<int>* electrons,
		std::vector<int> *muons, int* lepton_0, int* lepton_1, std::string mode) {
	/*Checks whether an OS lepton-lepton pair exists, and if so returns true and points lepton_0
	and lepton_1 to the selected particles*/
	std::vector<std::pair<int, int> > pairs; //Initialise array for OS pairs
	if (mode == "muons") {
		for (int m0 : *muons) { //Loop through muons
			for (int m1 : *muons) {
				if (m0 == m1) continue;
				if (reader->Muon_Charge[m0] != reader->Muon_Charge[m1]) {
					pairs.push_back(std::make_pair(m0, m1));
				}
			}
			if (pairs.size() == 1) { //Only one OS pair found
				*lepton_0 = pairs[0].first;
				*lepton_1 = pairs[0].second;
				if (reader->Muon_PT[*lepton_0] < reader->Muon_PT[*lepton_1]) { //Order muons by pT
					*lepton_1 = pairs[0].first;
					*lepton_0 = pairs[0].second;
				}
				return true;
			} else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
				double pT, highestPT = 0;
				std::pair<int, int> best;
				for (std::pair<int, int> p : pairs) { //Loop through pairs
					pT = std::abs(reader->Muon_PT[p.first])+std::abs(reader->Muon_PT[p.second]);
					if (pT > highestPT) {
						best = p;
						highestPT = pT;
					}
				}
				*lepton_0 = best.first;
				*lepton_1 = best.second;
				if (reader->Muon_PT[*lepton_0] < reader->Muon_PT[*lepton_1]) { //Order muons by pT
					*lepton_1 = best.first;
					*lepton_0 = best.second;
				}
				return true;
			} else { //No OS pairs found
				return false;
			}
		}
		return false;
	} else if (mode == "electrons") {
		for (int e0 : *electrons) { //Loop through electrons
			for (int e1 : *electrons) {
				if (e0 == e1) continue;
				if (reader->Electron_Charge[e0] != reader->Electron_Charge[e1]) {
					pairs.push_back(std::make_pair(e0, e1));
				}
			}
			if (pairs.size() == 1) { //Only one OS pair found
				*lepton_0 = pairs[0].first;
				*lepton_1 = pairs[0].second;
				if (reader->Electron_PT[*lepton_0] < reader->Electron_PT[*lepton_1]) { //Order electrons by pT
					*lepton_1 = pairs[0].first;
					*lepton_0 = pairs[0].second;
				}
				return true;
			} else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
				double pT, highestPT = 0;
				std::pair<int, int> best;
				for (std::pair<int, int> p : pairs) { //Loop through pairs
					pT = std::abs(reader->Electron_PT[p.first])+std::abs(reader->Electron_PT[p.second]);
					if (pT > highestPT) {
						best = p;
						highestPT = pT;
					}
				}
				*lepton_0 = best.first;
				*lepton_1 = best.second;
				if (reader->Electron_PT[*lepton_0] < reader->Electron_PT[*lepton_1]) { //Order electrons by pT
					*lepton_1 = best.first;
					*lepton_0 = best.second;
				}
				return true;
			} else { //No OS pairs found
				return false;
			}
		}
		return false;
	} else if (mode == "mixed") {
		for (int m : *muons) { //Loop through muons
			for (int e : *electrons) { //Loop through electrons
				if (reader->Muon_Charge[m] != reader->Electron_Charge[e]) {
					pairs.push_back(std::make_pair(m, e));
				}
			}
			if (pairs.size() == 1) { //Only one OS pair found
				*lepton_0 = pairs[0].first;
				*lepton_1 = pairs[0].second;
				return true;
			} else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
				double pT, highestPT = 0;
				std::pair<int, int> best;
				for (std::pair<int, int> p : pairs) { //Loop through pairs
					pT = std::abs(reader->Muon_PT[p.first])+std::abs(reader->Electron_PT[p.second]);
					if (pT > highestPT) {
						best = p;
						highestPT = pT;
					}
				}
				*lepton_0 = best.first;
				*lepton_1 = best.second;
				return true;
			} else { //No OS pairs found
				return false;
			}
		}
		return false;
	}
	return false;
}

bool getOSTauTauPair(delphesReader* reader, std::vector<int>* taus, int* tau_0, int* tau_1) {
	/*Checks whether an OS tau-tau pair exists, and if so returns true and points tau and lepton to
	the	selected particles*/
	std::vector<std::pair<int, int> > pairs; //Initialise array for OS tau pairs
	for (int t0 : *taus) { //Loop through taus
		for (int t1 : *taus) {
			if (t0 == t1) continue;
			if (reader->Jet_Charge[t0] != reader->Jet_Charge[t1]) {
				pairs.push_back(std::make_pair(t0, t1));
			}
		}
		if (pairs.size() == 1) { //Only one OS pair found
			*tau_0 = pairs[0].first;
			*tau_1 = pairs[0].second;
			if (reader->Jet_PT[*tau_0] < reader->Jet_PT[*tau_1]) { //Order taus by pT
				*tau_1 = pairs[0].first;
				*tau_0 = pairs[0].second;
			}
			return true;
		} else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
			double pT, highestPT = 0;
			std::pair<int, int> best;
			for (std::pair<int, int> p : pairs) { //Loop through pairs
				pT = std::abs(reader->Jet_PT[p.first])+std::abs(reader->Jet_PT[p.second]);
				if (pT > highestPT) {
					best = p;
					highestPT = pT;
				}
			}
			*tau_0 = best.first;
			*tau_1 = best.second;
			if (reader->Jet_PT[*tau_0] < reader->Jet_PT[*tau_1]) { //Order taus by pT
				*tau_1 = best.first;
				*tau_0 = best.second;
			}
			return true;
		} else { //No OS pairs found
			return false;
		}
	}
	return false;
}

int getNElectrons(delphesReader* reader) {
	/*Returns number of reco. electrons in event*/
	int nElectrons = 0;
	for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
		if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
				&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Quality electron
			nElectrons++;
		}
	}
	return nElectrons;
}

int getNMuons(delphesReader* reader) {
	/*Returns number of reco. muons in event*/
	int nMuons = 0;
	for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
		if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
				&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Quality muon
			nMuons++;
		}
	}
	return nMuons;
}

TLorentzVector getTauHadron(delphesReader* reader, int tau) {
	/*Returns 4-vector of hadronically decaying tau*/
	TLorentzVector t;
	t.SetPtEtaPhiM(reader->Jet_PT[tau], reader->Jet_Eta[tau], reader->Jet_Phi[tau], reader->Jet_Mass[tau]);
	return t;
}

TLorentzVector getTauLepton(delphesReader* reader, int lepton, std::string mode) {
	/*Returns 4-vector leptonically decaying tau*/
	TLorentzVector l;
	if (mode == "electron") {
		l.SetPtEtaPhiM(reader->Electron_PT[lepton], reader->Electron_Eta[lepton], reader->Electron_Phi[lepton], eMass);
	} else if (mode == "muon"){
		l.SetPtEtaPhiM(reader->Muon_PT[lepton], reader->Muon_Eta[lepton], reader->Muon_Phi[lepton], muMass);
	}
	return l;
}

TLorentzVector getBJet(delphesReader* reader, int bJet) {
	/*Returns 4-vector of b-jet*/
	TLorentzVector b;
	b.SetPtEtaPhiM(reader->Jet_PT[bJet], reader->Jet_Eta[bJet], reader->Jet_Phi[bJet], reader->Jet_Mass[bJet]);
	return b;
}

bool selectBJets(delphesReader* reader, std::vector<int>* bJets, int* bJet_0, int* bJet_1) {
	/*Checks is a pair of b-jets exists, returning true if so and pointing bJet_0 and bJet_1 to
	selected jets. Selects pair of jets invariant mass closest to 125 GeV*/
	if (bJets->size() == 2) { //Only two b jets found
		*bJet_0 = (*bJets)[0];
		*bJet_1 = (*bJets)[1];
		if (getBJet(reader, *bJet_0).Pt() < getBJet(reader, *bJet_1).Pt()) {
			*bJet_1 = (*bJets)[0];
			*bJet_0 = (*bJets)[1];
		}
		return true;
	} else if (bJets->size() > 2) { //More than two b jets: select pair with invariant mass closest to 125 GeV
		double deltaMin = -1;
		double delta;
		TLorentzVector jet_i, jet_j, jet_combined;
		int iMin, jMin;
		for (int i : *bJets) {
			jet_i = getBJet(reader, i);
			for (int j : *bJets) {
				if (i == j) continue;
				jet_j = getBJet(reader, j);
				jet_combined = jet_i + jet_j;
				delta = std::abs(125-jet_combined.M());
				if (deltaMin > delta || deltaMin < 0) {
					deltaMin = delta;
					iMin = i;
					jMin = j;
				}
			}
		}
		*bJet_0 = iMin;
		*bJet_1 = jMin;
		if (getBJet(reader, *bJet_0).Pt() < getBJet(reader, *bJet_1).Pt()) {
			*bJet_1 = iMin;
			*bJet_0 = jMin;
		}
		return true;
	} else { //Less than two b jets found
		return false;
	}
}

void makeDirs(std::string outputName) {
	/*Makes directory structure for saving outputs*/
	std::vector<std::string> dirs;
	dirs.push_back("../outputs");
	dirs.push_back("../outputs/" + outputName);
	for (std::string dir : dirs) {
		system((char*)("mkdir -p " + dir).c_str());
	}
}

std::vector<std::string> getInputs(std::string inputChar) {
	/*Filter input for file list*/
	std::string input(inputChar);
	const std::string dir = input.substr(0, input.rfind("/")+1);
	std::vector<std::string> inputFiles;
	const boost::regex mask(input.substr(input.rfind("/")+1));
	std::cout << "Getting input files from directory: " << dir << "\nUsing mask: " << mask << "\n";
	boost::filesystem::directory_iterator end;
	for (boost::filesystem::directory_iterator i(dir); i != end; ++i) {
		if (!boost::filesystem::is_regular_file(i->status())) continue;
		boost::smatch what;
		if(!boost::regex_match(i->path().filename().string(), what, mask)) continue;
		inputFiles.push_back(i->path().string());
	}
	return inputFiles;
}

TLorentzVector getHiggs2Taus(delphesReader* reader, TLorentzVector t_0, TLorentzVector t_1) {
	/*Returns 4-vector of Higgs->tau tau*/
	TLorentzVector higgs, mPT;
	mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0); //TODO Check this
	higgs = t_0 + t_1 + mPT;
	return higgs;
}

inline TLorentzVector getHiggs2Bs(TLorentzVector b_0, TLorentzVector b_1) {
	/*Returns 4-vector of Higgs->b b*/
	return b_0 + b_1;
}

inline TLorentzVector getDiHiggs(TLorentzVector higgs_0, TLorentzVector higgs_1) {
	/*Returns 4-vector of di_Higgs*/
	return higgs_0 + higgs_1;
}

int getBin(std::vector<double> edges, double value) {
	/*Returns bin number*/
	for (int i = 1; i < edges.size(); i++) {
		if (value < edges[i]) {
			return i-1;
		}
	}
	return edges.size()-1;
}

inline double getMT(double pT, double mPT, double dphi) {
	return sqrt(2*pT*mPT*(1-cos(dphi)));
}

void showHelp() {
	/*Show help for input arguments*/
	std::cout << "-i : input mask\n";
	std::cout << "-o : output name\n";
	std::cout << "-t : use MC truth cut [0/1], default 0\n";
	std::cout << "-s : run event selection [0/1], default 1\n";
	std::cout << "-d : run in debug mode [0/1], default 0\n";
	std::cout << "-m : output information for MVA selection [0/1], default 0\n";
}

void checkMother(int m, TClonesArray* particles, std::vector<int>* matches) {
	/*Checks to see whether particle came from tau*/
	GenParticle* mother = (GenParticle*)particles->At(m);
	if (std::abs(mother->PID) == 15) {
		matches->push_back(m); //If mother is a tau
	} else if (std::abs(mother->PID) == 24) { //If mother is a W boson
		if (mother->M1 > -1) {
			if (std::abs(((GenParticle*)particles->At(mother->M1))->PID) == 15) matches->push_back(mother->M1); //If grandmother is a tau
		}
		if (mother->M2 > -1) {
			if (std::abs(((GenParticle*)particles->At(mother->M2))->PID) == 15) matches->push_back(mother->M2); //If grandmother is a tau
		}
	}
}

std::vector<int> checkLepton(GenParticle* particle, TClonesArray* particles, int charge) {
	/*Checks to see whether lepton came from tau*/
	std::vector<int> matches;
	int tau = -1;
	if (!particle->IsPU && particle->Charge == charge) { //If particle is not PU and charge is correct
		if (particle->M1 > -1) checkMother(particle->M1, particles, &matches); //Check M1
		if (particle->M2 > -1) checkMother(particle->M2, particles, &matches); //Check M2
	}
	return matches;
}

std::vector<int> checkCommonMother(std::vector<int> a, std::vector<int> b, TClonesArray* particles) {
	/*Check to see if a pair of particle came from the decay of the same Higgs*/
	std::vector<int> pair;
	for (int i : a) {
		for (int j : b) {
			if (i == j) continue; //Same particle
			if (((GenParticle*)particles->At(i))->M1 != -1) {
				if (((GenParticle*)particles->At(i))->M1 == ((GenParticle*)particles->At(j))->M1 &&
						std::abs(((GenParticle*)particles->At(((GenParticle*)particles->At(i))->M1))->PID) == 25) {
					pair.push_back(i);
					pair.push_back(j);
					return pair;
				}
				if (((GenParticle*)particles->At(i))->M1 == ((GenParticle*)particles->At(j))->M2 &&
						std::abs(((GenParticle*)particles->At(((GenParticle*)particles->At(i))->M1))->PID) == 25) {
					pair.push_back(i);
					pair.push_back(j);
					return pair;
				}
			}
			if (((GenParticle*)particles->At(i))->M2 != -1) {
				if (((GenParticle*)particles->At(i))->M2 == ((GenParticle*)particles->At(j))->M1 &&
						std::abs(((GenParticle*)particles->At(((GenParticle*)particles->At(i))->M2))->PID) == 25) {
					pair.push_back(i);
					pair.push_back(j);
					return pair;
				}
				if (((GenParticle*)particles->At(i))->M2 == ((GenParticle*)particles->At(j))->M2 &&
						std::abs(((GenParticle*)particles->At(((GenParticle*)particles->At(i))->M2))->PID) == 25) {
					pair.push_back(i);
					pair.push_back(j);
					return pair;
				}
			}
		}
	}
	return pair;
}

std::string typeLookup(std::string mode) {
	/*Lookup for histogram bin name*/
	if (mode == "electron:electron") return "ee";
	if (mode == "muon:muon") return "#mu#mu";
	if (mode == "tau:tau") return "#tau_{h}#tau_{h}";
	if (mode == "electron:muon" || mode == "muon:electron") return "e#mu";
	if (mode == "electron:tau" || mode == "tau:electron") return "e#tau_{h}";
	if (mode == "tau:muon" || mode == "muon:tau") return "#mu#tau_{h}";
	return "error";
}

int moveToEnd(int p, TClonesArray* particles) {
	/*Selects particle at end of 'decay' chain in event record*/
	return p;
	GenParticle* mother = (GenParticle*)particles->At(p);
	while (((GenParticle*)particles->At(mother->D1))->PID == mother->PID &&
			((GenParticle*)particles->At(mother->D2))->PID == mother->PID) {
		p = mother->D1;
		mother = (GenParticle*)particles->At(p);
	}
	return p;
}

bool checkDiJet(TClonesArray* jets, TClonesArray* particles, int j_0, int j_1, int mother, int pID,
		int* swap, TH1D* dRPlot, double R) {
	/*Checks whether the particles are within their nearest jet*/
	//Associate particles to closest found jet___
	int p_0 = -1, p_1 = -1;
	Jet* jet_0 = (Jet*)jets->At(j_0);
	Jet* jet_1 = (Jet*)jets->At(j_1);
	GenParticle* higgs = (GenParticle*)particles->At(mother);
	if (std::abs(((GenParticle*)particles->At(moveToEnd(higgs->D1, particles)))->PID) != pID) { //Make sure decays products are correct
		std::cout << "Something's gone wrong in h->" + doubleToString(pID) + " -" + doubleToString(pID) + "!\n";
		return false;
	}
	if (((GenParticle*)particles->At(higgs->D1))->P4().DeltaR(jet_0->P4()) <
        ((GenParticle*)particles->At(higgs->D1))->P4().DeltaR(jet_1->P4())) {
		p_0 = higgs->D1;
        p_1 = higgs->D2;
        if (swap != NULL) *swap = 0;
    } else {
        p_0 = higgs->D2;
		p_1 = higgs->D1;
        if (swap != NULL) *swap = 1;
    }
	//___________________________________________
	//Check jets_________________________________
	double dR_0 = ((GenParticle*)particles->At(p_0))->P4().DeltaR(jet_0->P4());
	double dR_1 = ((GenParticle*)particles->At(p_1))->P4().DeltaR(jet_1->P4());
	if (dR_0 > R || dR_1 > R) { //particle(s) outside jet
		return false;
	}
	//___________________________________________
	//Accept association and fill plot___________
	dRPlot->Fill(dR_0);
	dRPlot->Fill(dR_1);
	//___________________________________________
	return true;
}

int ancestrySearch(GenParticle* child, GenParticle* parent_0, GenParticle* parent_1, TClonesArray* particles) {
	/*Recursive search through child's ancestry for parent 0 or 1. If found returns 0 or 1. If not found returns -1*/
	int ancestor = -1;
	GenParticle* mother;
	if (child->M1 > 0) { //Check first mother
		mother = (GenParticle*)particles->At(child->M1);
		if (mother->GetUniqueID() == parent_0->GetUniqueID()) {
			return 0;
		} else if (mother->GetUniqueID() == parent_1->GetUniqueID()) {
			return 1;
		} else {
			ancestor = ancestrySearch(mother, parent_0, parent_1, particles); //Recursive search
		}
	}
	if (ancestor == -1 && child->M2 > 0) { //Check second mother
		mother = (GenParticle*)particles->At(child->M2);
		if (mother->GetUniqueID() == parent_0->GetUniqueID()) {
			return 0;
		} else if (mother->GetUniqueID() == parent_1->GetUniqueID()) {
			return 1;
		} else {
			ancestor = ancestrySearch(mother, parent_0, parent_1, particles); //Recursive search
		}
	}
	return ancestor;
}

bool correctDecayChannel(std::string input, Long64_t cEvent,
		std::map<std::string, TH1D*>* plots=NULL, int* hBB=NULL, int* hTauTau=NULL) {
	/*Make sure event is hh->bbtautau, and point hbb and htautau to the Higgs*/
	TChain *chain = new TChain("Delphes");
	chain->Add(input.c_str());
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	treeReader->ReadEntry(cEvent);
	bool hBBFound = false, hTauTauFound = false;
	int nHiggs = 0;
	if (plots != NULL) (*plots)["cuts"]->Fill("hh->bb#tau#tau check", 1);
	for (int p = 0; p < branchParticle->GetEntriesFast(); ++p) {
		if (std::abs(((GenParticle*)branchParticle->At(p))->PID) == 25) { //Particle is Higgs
			if (((GenParticle*)branchParticle->At(p))->D1 >= 0 && ((GenParticle*)branchParticle->At(p))->D2 >= 0) { //Daughters exists
				if (((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID != 25 &&
						((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID != 25) {
					nHiggs++;
					if (plots != NULL) (*plots)["higgsDecay"]->Fill(std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID));
					if (plots != NULL) (*plots)["higgsDecay"]->Fill(std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID));
					if (std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID) == 5
							&& std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID) == 5) { //Daughters are b quarks
						hBBFound = true;
						if (hBB != NULL) *hBB = p; //Point to Higgs
						if (hBBFound && hTauTauFound) { //h->bb and h->tautau found, so accept event
							if (plots != NULL) (*plots)["cuts"]->Fill("hh->bb#tau#tau pass", 1);
							chain->Delete();
							delete treeReader;
							return true;
						}
					}
					if (std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID) == 15
							&& std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID) == 15) { //Daughters are taus
						hTauTauFound = true;
						if (hTauTau != NULL) *hTauTau = p; //Point to Higgs
						if (hBBFound && hTauTauFound) { //h->bb and h->tautau found, so accept event
							if (plots != NULL) (*plots)["cuts"]->Fill("hh->bb#tau#tau pass", 1);
							chain->Delete();
							delete treeReader;
							return true;
						}
					}
				}
			}
			if (nHiggs >= 2) break; //Both Higgs found
		}
	}
	chain->Delete();
	delete treeReader;
	return false; //Both h->bb and h->tautau not found
}

std::map<std::string, std::string> getOptions(int argc, char* argv[]) {
	/*Interpret input arguments*/
	std::map<std::string, std::string> options;
	options.insert(std::make_pair("-i", "")); //Input mask
	options.insert(std::make_pair("-o", "")); //Output name
	options.insert(std::make_pair("-t", "0")); //MC truth cut
	options.insert(std::make_pair("-s", "1")); //Event selection
	options.insert(std::make_pair("-d", "0")); //Debug mode
	options.insert(std::make_pair("-m", "0")); //Output information for MVA selection
	if (argc >= 2) {
		std::string option(argv[1]);
		if (option == "-h" || option == "--help") {
			showHelp();
			options.clear();
			return options;
		}
	}
	for (int i = 1; i < argc; i = i+2) {
		std::string option(argv[i]);
		std::string argument(argv[i+1]);
		if (option == "-h" || option == "--help" || argument == "-h" || argument == "--help") {
			showHelp();
			options.clear();
			return options;
		}
		options[option] = argument;
	}
	if (options["-i"] == "" || options["-o"] == "") {
		showHelp();
		options.clear();
		return options;
	}
	if (options["-d"] == "1") {
		debug = true;
		std::cout << "Running in debug mode\n";
	}
	return options;
}

TMatrixD decomposeVector(TLorentzVector in) {
	TMatrixD out(3, 3);
	out(0, 0) = in.Px()*in.Px();
	out(0, 1) = in.Px()*in.Py();
	out(0, 2) = in.Px()*in.Pz();
	out(1, 0) = in.Py()*in.Px();
	out(1, 1) = in.Py()*in.Py();
	out(1, 2) = in.Py()*in.Pz();
	out(2, 0) = in.Pz()*in.Px();
	out(2, 1) = in.Pz()*in.Py();
	out(2, 2) = in.Pz()*in.Pz();
	return out;
}

void appendSphericity(TMatrixD* mat, double* div, TLorentzVector mom) {
	TMatrixD decomp = decomposeVector(mom);
	*mat += decomp;
	*div += pow(mom.P(), 2);
}

void appendSpherocity(TMatrixD* mat, double* div, TLorentzVector mom) {
	TMatrixD decomp = decomposeVector(mom);
	decomp *= 1/std::abs(mom.P());
	*mat += decomp;
	*div += std::abs(mom.P());
}

std::vector<double> getEigenValues(TMatrixD in) {
	/*Return vector of sorted, nomalised eigenvalues of parssed matrix*/
	TMatrixD eigenMatrix = TMatrixDEigen(in).GetEigenValues();
	std::vector<double> eigenValues(3);
	eigenValues[0] = eigenMatrix(0, 0);
	eigenValues[1] = eigenMatrix(1, 1);
	eigenValues[2] = eigenMatrix(2, 2);
	std::sort(eigenValues.begin(), eigenValues.end(), std::greater<double>());
	double sum = 0;
	for (double n : eigenValues)
		sum += n;
	std::for_each(eigenValues.begin(), eigenValues.end(), [sum](double i) { return i/sum; });
	return eigenValues;
}

void getEventShapes(std::vector<double> sphericityV, std::vector<double> spherocityV,
		double* sphericity, double* spherocity,
		double* aplanarity, double* aplanority,
		double* upsilon, double* dShape) {
	*sphericity = (3/2)*(sphericityV[1]+sphericityV[2]);
	*spherocity = (3/2)*(spherocityV[1]+spherocityV[2]);
	*aplanarity = 3*sphericityV[2]/2;
	*aplanority = 3*spherocityV[2]/2;
	*upsilon = sqrt(3.0)*(sphericityV[1]-sphericityV[2])/2;
	*dShape = 27*spherocityV[0]*spherocityV[1]*spherocityV[2];
}

void getGlobalEventInfo(std::string input, Long64_t cEvent,
		double*  hT, double*  sT, double* centrality, double* eVis,
		int* nJets, int* nBJets, int* nTauJets, int* nPhotons,
		double* minJetPT, double* meanJetPT, double* maxJetPT,
		double* minJetMass, double* meanJetMass, double* maxJetMass,
		double* minJetEta, double* meanJetEta, double* maxJetEta,
		double* sphericity, double* spherocity,
		double* aplanarity, double* aplanority,
		double* upsilon, double* dShape) {
	/*Fills referenced variables with global event information*/
	if (debug) std::cout << "Getting global event info\n";
	//Load event info____________________________
	TChain *chain = new TChain("Delphes");
	chain->Add(input.c_str());
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
	treeReader->ReadEntry(cEvent);
	if (debug) std::cout << "Loaded info\n";
	//___________________________________________
	//Reset variables____________________________
	*hT = 0;
	*sT = 0;
	*centrality = 0;
	*eVis = 0;
	*nJets = 0;
	*nBJets = 0;
	*nTauJets = 0;
	*nPhotons = 0;
	*minJetPT = -1;
	*meanJetPT = 0;
	*maxJetPT = -1;
	*minJetMass = -1;
	*meanJetMass = 0;
	*maxJetMass = -1;
	*minJetEta = -1;
	*meanJetEta = 0;
	*maxJetEta = -1;
	*sphericity = 0;
	*spherocity = 0;
	*aplanarity = 0;
	*aplanority = 0;
	*upsilon = 0;
	*dShape = 0;
	//___________________________________________
	//Initialise holders_________________________
	Electron* electron;
	Muon* muon;
	Jet* jet;
	MissingET* mPT;
	Photon* photon;
	TMatrixD sphericityT(3, 3), spherocityT(3, 3);
	double sphericityD = 0, spherocityD = 0;
	//___________________________________________
	//Loop through objects_____________________
	for (int i = 0; i < branchElectron->GetEntriesFast(); ++i) { //Loop over all electrons in event
		electron = (Electron*)branchElectron->At(i);
		*sT += electron->PT;
		*eVis += electron->P4().E();
		*centrality += electron->PT;
		appendSphericity(&sphericityT, &sphericityD, electron->P4());
		appendSpherocity(&spherocityT, &spherocityD, electron->P4());
	}
	for (int i = 0; i < branchPhoton->GetEntriesFast(); ++i) { //Loop over all photons in event
		photon = (Photon*)branchPhoton->At(i);
		if (photon->Particles.GetEntriesFast() != 1) continue; //Skip photons with references to multiple particles
		*nPhotons += 1;
		*sT += photon->PT;
		*eVis += photon->P4().E();
		*centrality += photon->PT;
		appendSphericity(&sphericityT, &sphericityD, photon->P4());
		appendSpherocity(&spherocityT, &spherocityD, photon->P4());
	}
	for (int i = 0; i < branchMuon->GetEntriesFast(); ++i) { //Loop over all muons in event
		muon = (Muon*)branchMuon->At(i);
		*sT += muon->PT;
		*eVis += muon->P4().E();
		*centrality += muon->PT;
		appendSphericity(&sphericityT, &sphericityD, muon->P4());
		appendSpherocity(&spherocityT, &spherocityD, muon->P4());
	}
	mPT = (MissingET*)branchMissingET->At(0);
	*sT += mPT->MET;
	for (int i = 0; i < branchJet->GetEntriesFast(); ++i) { //Loop over all jets in event
		jet = (Jet*)branchJet->At(i);
		*hT += sqrt(pow(jet->Mass, 2)+pow(jet->PT, 2));
		*sT += sqrt(pow(jet->Mass, 2)+pow(jet->PT, 2));
		*centrality += jet->PT;
		*eVis += jet->P4().E();
		*nJets += 1;
		if (jet->TauTag) *nTauJets += 1;
		if (jet->BTag) *nBJets += 1;
		if (*minJetPT == -1 | jet->PT < *minJetPT) *minJetPT = jet->PT;
		*meanJetPT += jet->PT;
		if (jet->PT > *maxJetPT) *maxJetPT = jet->PT;
		if (*minJetMass == -1 | jet->Mass < *minJetMass) *minJetMass = jet->Mass;
		*meanJetMass += jet->Mass;
		if (jet->Mass > *maxJetMass) *maxJetMass = jet->Mass;
		if (*minJetEta == -1 | std::abs(jet->Eta) < *minJetEta) *minJetEta = std::abs(jet->Eta);
		*meanJetEta += jet->Eta;
		if (std::abs(jet->Eta) > *maxJetEta) *maxJetEta = std::abs(jet->Eta);
		appendSphericity(&sphericityT, &sphericityD, jet->P4());
		appendSpherocity(&spherocityT, &spherocityD, jet->P4());
	}
	//___________________________________________
	//Finalise variabales________________________
	*centrality /= *eVis;
	*meanJetPT /= *nJets;
	*meanJetMass /= *nJets;
	*meanJetEta /= *nJets;
	sphericityT *= 1/sphericityD;
	spherocityT *= 1/spherocityD;
	//___________________________________________
	//Calculate event shapes_____________________
	if (debug) std::cout << "Calculating global event shapes\n";
	std::vector<double> sphericityV = getEigenValues(sphericityT);
	std::vector<double> spherocityV = getEigenValues(spherocityT);
	getEventShapes(sphericityV, spherocityV,
			sphericity, spherocity,
			aplanarity, aplanority,
			upsilon, dShape);
	//___________________________________________
	chain->Delete();
	delete treeReader;
}

void getPrimaryEventShapes(TLorentzVector v_tau_0, TLorentzVector v_tau_1, TLorentzVector v_bJet_0, TLorentzVector v_bJet_1,
		double* sphericity, double* spherocity,
		double* aplanarity, double* aplanority,
		double* upsilon, double* dShape) {
	/*Sets values of referenced event-shape variables for final-states*/
	if (debug) std::cout << "Getting primary event shapes\n";
	TMatrixD sphericityT(3, 3), spherocityT(3, 3);
	double sphericityD = 0, spherocityD = 0;
	//Populate tensors___________________________
	appendSphericity(&sphericityT, &sphericityD, v_tau_0);
	appendSpherocity(&spherocityT, &spherocityD, v_tau_0);
	appendSphericity(&sphericityT, &sphericityD, v_tau_1);
	appendSpherocity(&spherocityT, &spherocityD, v_tau_1);
	appendSphericity(&sphericityT, &sphericityD, v_bJet_0);
	appendSpherocity(&spherocityT, &spherocityD, v_bJet_0);
	appendSphericity(&sphericityT, &sphericityD, v_bJet_1);
	appendSpherocity(&spherocityT, &spherocityD, v_bJet_1);
	sphericityT *= 1/sphericityD;
	spherocityT *= 1/spherocityD;
	//___________________________________________
	//Calculate event shapes_____________________
	if (debug) std::cout << "Calculating primary event shapes\n";
	std::vector<double> sphericityV = getEigenValues(sphericityT);
	std::vector<double> spherocityV = getEigenValues(spherocityT);
	getEventShapes(sphericityV, spherocityV,
			sphericity, spherocity,
			aplanarity, aplanority,
			upsilon, dShape);
	//___________________________________________
}

bool truthCut(std::string input, Long64_t cEvent, int b_0, int b_1, int l_0, int l_1, int hBB, int hTauTau,
		std::string mode, std::map<std::string, TH1D*>* plots,
		TLorentzVector* v_gen_higgs_bb, TLorentzVector* v_gen_higgs_tt,
		TLorentzVector* v_gen_tau_0, TLorentzVector* v_gen_tau_1,
		TLorentzVector* v_gen_bJet_0, TLorentzVector* v_gen_bJet_1) {
	/*Checks whether selected final states are correct*/
	double jetRadius = 0.5;
	TChain *chain = new TChain("Delphes");
	chain->Add(input.c_str());
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	treeReader->ReadEntry(cEvent);
	if (debug) std::cout << "Loading data for MC truth cut on event mode " << mode << "\n";
	//Check if selected final states are correct_
	(*plots)["cuts"]->Fill("MC-truth check", 1);
	int swap;
	//Check b jets_______________________________
	GenParticle *bJet_0, *bJet_1;
	GenParticle* higgs = (GenParticle*)branchParticle->At(hBB);
	(*plots)["cuts"]->Fill("b-jets check", 1);
	if (debug) std::cout << "Checking b-jets\n";
	if (!checkDiJet(branchJet, branchParticle, b_0, b_1, hBB, 5, &swap, (*plots)["bMatch"], jetRadius)) {
		if (debug) std::cout << "MC check fails due to di-Jet on b-jets check\n";
		chain->Delete();
		return false; //b-jet selection incorrect
	}
	if (debug) std::cout << "Both b jets confirmed\n";
	(*plots)["cuts"]->Fill("b-jets pass", 1);
	if (swap) {
		bJet_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		bJet_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
	} else {
		bJet_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		bJet_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
	}
	//___________________________________________
	//Check taus_________________________________
	if (debug) std::cout << "Checking taus\n";
	std::vector<std::string> options;
	boost::split(options, mode, boost::is_any_of(":"));
	(*plots)["cuts"]->Fill("#taus check", 1);
	(*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(mode) + " check").c_str(), 1);
	GenParticle *tau_0, *tau_1;
	higgs = (GenParticle*)branchParticle->At(hTauTau);
	if (options[0] == "tau" && options[1] == "tau") {
		//h->tau_h tau_h_________________________
		if (!checkDiJet(branchJet, branchParticle, l_0, l_1, hTauTau, 15, &swap, (*plots)["tauMatch"], jetRadius)) {
			if (debug) std::cout << "MC check fails due to di-Jet on tau-jets check\n";
			chain->Delete();
			return false; //tau-jet selection incorrect
		}
		if (swap) {
			tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
			tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		} else {
			tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
			tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		}
		(*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(mode) + " pass").c_str(), 1);
		//_______________________________________
	} else if ((options[0] == "tau" && options[1] != "tau") || (options[0] != "tau" && options[1] == "tau")) {
		//h->tau_h light-lepton__________________
		//Load objects___________________________
		tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		GenParticle* lightLepton;
		Jet* tauJet;
		for (int i = 0; i < 2; i++) {
			int l = l_0;
			if (i == 1) {
				l = l_1;
			}
			if (options[i] == "tau") {
				tauJet = (Jet*)branchJet->At(l);
			} else if (options[i] == "muon") {
				lightLepton = (GenParticle*)((Muon*)branchMuon->At(l))->Particle.GetObject();
			} else if (options[i] == "electron") {
				lightLepton = (GenParticle*)((Electron*)branchElectron->At(l))->Particle.GetObject();
			}
		}
		//_______________________________________
		//Check taus_____________________________
		int leptonMother = ancestrySearch(lightLepton, tau_0, tau_1, branchParticle);
		if (leptonMother == -1) {
			if (debug) std::cout << "MC check fails due to ancestry check\n";
			chain->Delete();
			return false; //Light lepton did not come from tau decay
		}
		GenParticle* tauh;
		if (leptonMother == 0) {
			tauh = tau_1;
			tau_1 = tau_0; //Reassociate 0 to tau and 1 to lepton
			tau_0 = tauh;
		} else {
			tauh = tau_0;
		}
		if (tauh->P4().DeltaR(tauJet->P4()) > jetRadius) {
			if (debug) std::cout << "MC check fails due to tau-jet check\n";
			chain->Delete();
			return false; //Tau outside selected jet
		}
		(*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(mode) + " pass").c_str(), 1);
		//_______________________________________
		//_______________________________________
	} else {
		//h->light-lepton light-lepton___________
		//Load objects___________________________
		GenParticle* higgs = (GenParticle*)branchParticle->At(hTauTau);
		tau_0 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D2, branchParticle));
		GenParticle *lightLepton_0, *lightLepton_1;
		if (options[0] == "muon") {
			lightLepton_0 = (GenParticle*)((Muon*)branchMuon->At(l_0))->Particle.GetObject();
		} else if (options[0] == "electron") {
			lightLepton_0 = (GenParticle*)((Electron*)branchElectron->At(l_0))->Particle.GetObject();
		}
		if (options[1] == "muon") {
			lightLepton_1 = (GenParticle*)((Muon*)branchMuon->At(l_1))->Particle.GetObject();
		} else if (options[1] == "electron") {
			lightLepton_1 = (GenParticle*)((Electron*)branchElectron->At(l_1))->Particle.GetObject();
		}
		//_______________________________________
		//Check taus_____________________________
		int leptonMother_0 = ancestrySearch(lightLepton_0, tau_0, tau_1, branchParticle);
		if (leptonMother_0 == -1) {
			if (debug) std::cout << "MC check fails due to ancestry check\n";
			chain->Delete();
			return false; //Light lepton 0 did not come from tau decay
		}
		int leptonMother_1 = ancestrySearch(lightLepton_1, tau_0, tau_1, branchParticle);
		if (leptonMother_1 == -1) {
			if (debug) std::cout << "MC check fails due to ancestry check\n";
			chain->Delete();
			return false; //Light lepton 1 did not come from tau decay
		}
		if (leptonMother_0 == leptonMother_1) {
			if (debug) std::cout << "MC check fails due to both leptons coming from same tau\n";
			chain->Delete();
			return false; //Leptons both came from same mother (somehow)
		}
		(*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(mode) + " pass").c_str(), 1);
		if ((lightLepton_0->PT > lightLepton_1->PT & leptonMother_0 == 1) |
				(lightLepton_0->PT < lightLepton_1->PT & leptonMother_0 == 0)) {
			tau_0 = tau_1;
			tau_1 = (GenParticle*)branchParticle->At(moveToEnd(higgs->D1, branchParticle));
		}
		//_______________________________________
		//_______________________________________
	}
	if (debug) std::cout << "Both taus confirmed\n";
	//___________________________________________
	if (debug) std::cout << "Event accepted\n";
	(*plots)["cuts"]->Fill("#taus pass", 1);
	(*plots)["cuts"]->Fill("MC-truth pass", 1);
	//Get vectors for regression_________________
	*v_gen_higgs_bb = ((GenParticle*)branchParticle->At(hBB))->P4();
	*v_gen_higgs_tt = ((GenParticle*)branchParticle->At(hTauTau))->P4();
	*v_gen_tau_0 = tau_0->P4();
	*v_gen_tau_1 = tau_1->P4();
	*v_gen_bJet_0 = bJet_0->P4();
	*v_gen_bJet_1 = bJet_1->P4();
	//___________________________________________
	chain->Delete();
	delete treeReader;
	return true;
}

void printMVASelectionInfo(std::string input, std::string output, std::string truthCut) {
	/*Output 4-vectors of reconstructed object for an MVA-based selection*/
	//Initialise variables_______________________
	std::cout << "Initialising variables\n";
	std::vector<std::string> inputs = getInputs(input);
	std::cout << inputs.size() << " input files found\n";
	Electron *electron;
	Photon *photon;
	Muon *muon;
	Jet *jet;
	MissingET *met;
	int hBB, hTauTau;
	double l_0_eta, l_1_eta, l_2_eta;
	double l_0_phi, l_1_phi, l_2_phi;
	double l_0_pt, l_1_pt, l_2_pt;
	double l_0_mass, l_1_mass, l_2_mass;
	double l_0_charge, l_1_charge, l_2_charge;
	int l_0_muon, l_1_muon, l_2_muon;
	double j_0_eta, j_1_eta, j_2_eta, j_3_eta, j_4_eta;
	double j_0_phi, j_1_phi, j_2_phi, j_3_phi, j_4_phi;
	double j_0_pt, j_1_pt, j_2_pt, j_3_pt, j_4_pt;
	double j_0_mass, j_1_mass, j_2_mass, j_3_mass, j_4_mass;
	double j_0_bTag, j_1_bTag, j_2_bTag, j_3_bTag, j_4_bTag;
	double j_0_tauTag, j_1_tauTag, j_2_tauTag, j_3_tauTag, j_4_tauTag;
	double met_phi, met_pt;
	double p_eta,  p_phi, p_pt;
	double ht_all, ht_l, ht_j;
	int nLeptons, nJets, nPhotons;
	double gen_hh_mass, gen_hBB_mass, gen_hTauTau_mass;
	TTree* data = new TTree("data", "data");
	data->Branch("l_0_eta", &l_0_eta);
	data->Branch("l_1_eta", &l_1_eta);
	data->Branch("l_2_eta", &l_2_eta);
	data->Branch("l_0_phi", &l_0_eta);
	data->Branch("l_1_phi", &l_1_eta);
	data->Branch("l_2_phi", &l_2_eta);
	data->Branch("l_0_pt", &l_0_pt);
	data->Branch("l_1_pt", &l_1_pt);
	data->Branch("l_2_pt", &l_2_pt);
	data->Branch("l_0_mass", &l_0_mass);
	data->Branch("l_1_mass", &l_1_mass);
	data->Branch("l_2_mass", &l_2_mass);
	data->Branch("l_0_charge", &l_0_charge);
	data->Branch("l_1_charge", &l_1_charge);
	data->Branch("l_2_charge", &l_2_charge);
	data->Branch("l_0_muon", &l_0_muon);
	data->Branch("l_1_muon", &l_1_muon);
	data->Branch("l_2_muon", &l_2_muon);
	data->Branch("j_0_eta", &j_0_eta);
	data->Branch("j_1_eta", &j_1_eta);
	data->Branch("j_2_eta", &j_2_eta);
	data->Branch("j_3_eta", &j_3_eta);
	data->Branch("j_4_eta", &j_4_eta);
	data->Branch("j_0_phi", &j_0_phi);
	data->Branch("j_1_phi", &j_1_phi);
	data->Branch("j_2_phi", &j_2_phi);
	data->Branch("j_3_phi", &j_3_phi);
	data->Branch("j_4_phi", &j_4_phi);
	data->Branch("j_0_pt", &j_0_pt);
	data->Branch("j_1_pt", &j_1_pt);
	data->Branch("j_2_pt", &j_2_pt);
	data->Branch("j_3_pt", &j_3_pt);
	data->Branch("j_4_pt", &j_4_pt);
	data->Branch("j_0_mass", &j_0_mass);
	data->Branch("j_1_mass", &j_1_mass);
	data->Branch("j_2_mass", &j_2_mass);
	data->Branch("j_3_mass", &j_3_mass);
	data->Branch("j_4_mass", &j_4_mass);
	data->Branch("j_0_bTag", &j_0_bTag);
	data->Branch("j_1_bTag", &j_1_bTag);
	data->Branch("j_2_bTag", &j_2_bTag);
	data->Branch("j_3_bTag", &j_3_bTag);
	data->Branch("j_4_bTag", &j_4_bTag);
	data->Branch("j_0_tauTag", &j_0_tauTag);
	data->Branch("j_1_tauTag", &j_1_tauTag);
	data->Branch("j_2_tauTag", &j_2_tauTag);
	data->Branch("j_3_tauTag", &j_3_tauTag);
	data->Branch("j_4_tauTag", &j_4_tauTag);
	data->Branch("met_phi", &met_phi);
	data->Branch("met_pt", &met_pt);
	data->Branch("p_eta", &p_eta);
	data->Branch("p_phi", &p_phi);
	data->Branch("p_pt", &p_pt);
	data->Branch("ht_all", &ht_all);
	data->Branch("ht_j", &ht_j);
	data->Branch("ht_l", &ht_l);
	data->Branch("nLeptons", &nLeptons);
	data->Branch("nJets", &nJets);
	data->Branch("nPhotons", &nPhotons);
	data->Branch("gen_hh_mass", &gen_hh_mass);
	data->Branch("gen_hBB_mass", &gen_hBB_mass);
	data->Branch("gen_hTauTau_mass", &gen_hTauTau_mass);
	//_______________________________________
	//Loop through input files_______________
	for (int f = 0; f < inputs.size(); f++) {
		//Load data__________________________
		std::cout << "Loading file " << f+1 << " of " << inputs.size() << "\n";
		TChain *chain = new TChain("Delphes");
		chain->Add(inputs[f].c_str());
		ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
		TClonesArray *branchParticle = treeReader->UseBranch("Particle");
		TClonesArray *branchElectron = treeReader->UseBranch("Electron");
		TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
		TClonesArray *branchMuon = treeReader->UseBranch("Muon");
		TClonesArray *branchJet = treeReader->UseBranch("Jet");
		TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
		//___________________________________
		//Loop through events________________
		Long64_t nEvents = treeReader->GetEntries();
		for (Long64_t entry = 0; entry < nEvents; ++entry) { //Loop over all events
			treeReader->ReadEntry(entry); //Load selected branches with data from specified event
			if (entry % 1000 == 0) std::cout << "Loop: " << entry << "/" << nEvents << ", " << 100*entry/nEvents << "%\n";
			if (truthCut == "1") {
				if (!correctDecayChannel(inputs[f], entry, NULL, &hBB, &hTauTau)) continue; //Event is not h->bbtautau
			}
			//Load objects and sort by pT________
			std::vector<std::pair<TLorentzVector, std::vector<double> > > leptons;
			std::vector<std::pair<TLorentzVector, double> > photons(1);
			std::vector<std::pair<Jet*, double> > jets(5);
			nLeptons = 0;
			nPhotons = 0;
			nJets = 0;
			ht_all = 0;
			ht_j = 0;
			ht_l = 0;
			for (int i = 0; i < branchElectron->GetEntriesFast(); ++i) { //Loop over all electrons in event
				electron = (Electron*)branchElectron->At(i);
				leptons.push_back(std::make_pair(electron->P4(), std::vector<double>{(double)electron->PT, (double)electron->Charge, -1.0, (double)electron->IsolationVar}));
				nLeptons++;
				ht_l += electron->PT;
				ht_all += electron->PT;
			}
			for (int i = 0; i < branchPhoton->GetEntriesFast(); ++i) { //Loop over all photons in event
				photon = (Photon*)branchPhoton->At(i);
				if (photon->Particles.GetEntriesFast() != 1) continue; //Skip photons with references to multiple particles
				photons.push_back(std::make_pair(photon->P4(), photon->PT));
				nPhotons++;
				ht_all += photon->PT;
			}
			for (int i = 0; i < branchMuon->GetEntriesFast(); ++i) { //Loop over all muons in event
				muon = (Muon*)branchMuon->At(i);
				leptons.push_back(std::make_pair(muon->P4(), std::vector<double>{(double)muon->PT, (double)muon->Charge, 1.0, (double)muon->IsolationVar}));
				nLeptons++;
				ht_l += muon->PT;
				ht_all += muon->PT;
			}
			for (int i = 0; i < branchJet->GetEntriesFast(); ++i) { //Loop over all jets in event
				jet = (Jet*)branchJet->At(i);
				jets.push_back(std::make_pair(jet, jet->PT));
				nJets++;
				ht_j += jet->PT;
				ht_all += jet->PT;
			}
			met = (MissingET*)branchMissingET->At(0);
			while (leptons.size() < 3) {
				leptons.push_back(std::make_pair(TLorentzVector(), std::vector<double>(1)));
			}
			std::sort(leptons.begin(), leptons.end(),
					[](const std::pair<TLorentzVector, std::vector<double> >& a, const std::pair<TLorentzVector, std::vector<double> >& b)
					{return a.second[0] > b.second[0];});
			std::sort(photons.begin(), photons.end(),
					[](const std::pair<TLorentzVector, double>& a, const std::pair<TLorentzVector, double>& b)
					{return a.second > b.second;});
			std::sort(jets.begin(), jets.end(),
					[](const std::pair<Jet*, double>& a, const std::pair<Jet*, double>& b)
					{return a.second > b.second;});
			//___________________________________
			//Select objects and fill branch_____
			l_0_eta = leptons[0].first.Eta();
			l_1_eta = leptons[1].first.Eta();
			l_2_eta = leptons[2].first.Eta();
			l_0_phi = leptons[0].first.Phi();
			l_1_phi = leptons[1].first.Phi();
			l_2_phi = leptons[2].first.Phi();
			l_0_pt = leptons[0].first.Pt();
			l_1_pt = leptons[1].first.Pt();
			l_2_pt = leptons[2].first.Pt();
			l_0_mass = (leptons[0].second[2] == 1.0)? muMass : (leptons[0].second[2] == -1.0)? eMass : 0.0;
			l_1_mass = (leptons[1].second[2] == 1.0)? muMass : (leptons[1].second[2] == -1.0)? eMass : 0.0;
			l_2_mass = (leptons[2].second[2] == 1.0)? muMass : (leptons[2].second[2] == -1.0)? eMass : 0.0;
			l_0_muon = leptons[0].second[2];
			l_1_muon = leptons[1].second[2];
			l_2_muon = leptons[2].second[2];
			l_0_charge = leptons[0].second[1];
			l_1_charge = leptons[1].second[1];
			l_2_charge = leptons[2].second[1];
			j_0_eta = (jets[0].second > 0)? jets[0].first->Eta : 0;
			j_1_eta = (jets[1].second > 0)? jets[1].first->Eta : 0;
			j_2_eta = (jets[2].second > 0)? jets[2].first->Eta : 0;
			j_3_eta = (jets[3].second > 0)? jets[3].first->Eta : 0;
			j_4_eta = (jets[4].second > 0)? jets[4].first->Eta : 0;
			j_0_phi = (jets[0].second > 0)? jets[0].first->Phi : 0;
			j_1_phi = (jets[0].second > 0)? jets[0].first->Phi : 0;
			j_2_phi = (jets[1].second > 0)? jets[1].first->Phi : 0;
			j_3_phi = (jets[2].second > 0)? jets[2].first->Phi : 0;
			j_4_phi = (jets[3].second > 0)? jets[3].first->Phi : 0;
			j_0_pt = jets[0].second;
			j_1_pt = jets[1].second;
			j_2_pt = jets[2].second;
			j_3_pt = jets[3].second;
			j_4_pt = jets[4].second;
			j_0_mass = (jets[0].second > 0)? jets[0].first->Mass : 0;
			j_1_mass = (jets[1].second > 0)? jets[1].first->Mass : 0;
			j_2_mass = (jets[2].second > 0)? jets[2].first->Mass : 0;
			j_3_mass = (jets[3].second > 0)? jets[3].first->Mass : 0;
			j_4_mass = (jets[4].second > 0)? jets[4].first->Mass : 0;
			j_0_bTag = (jets[0].second > 0)? jets[0].first->BTag : 0;
			j_1_bTag = (jets[1].second > 0)? jets[1].first->BTag : 0;
			j_2_bTag = (jets[2].second > 0)? jets[2].first->BTag : 0;
			j_3_bTag = (jets[3].second > 0)? jets[3].first->BTag : 0;
			j_4_bTag = (jets[4].second > 0)? jets[4].first->BTag : 0;
			j_0_tauTag = (jets[0].second > 0)? jets[0].first->TauTag : 0;
			j_1_tauTag = (jets[1].second > 0)? jets[1].first->TauTag : 0;
			j_2_tauTag = (jets[2].second > 0)? jets[2].first->TauTag : 0;
			j_3_tauTag = (jets[3].second > 0)? jets[3].first->TauTag : 0;
			j_4_tauTag = (jets[4].second > 0)? jets[4].first->TauTag : 0;
			met_phi = met->Phi;
			met_pt = met->MET;
			p_eta = (photons[0].second > 0)? photons[0].first.Eta() : 0;
			p_phi = (photons[0].second > 0)? photons[0].first.Phi() : 0;
			p_pt = (photons[0].second > 0)? photons[0].first.Pt() : 0;
			gen_hh_mass = getDiHiggs(((GenParticle*)branchParticle->At(hBB))->P4(),
					((GenParticle*)branchParticle->At(hTauTau))->P4()).M();
			gen_hBB_mass = ((GenParticle*)branchParticle->At(hBB))->P4().M();
			gen_hTauTau_mass = ((GenParticle*)branchParticle->At(hTauTau))->P4().M();
			data->Fill();
			//___________________________________
		}
		//_______________________________________
		chain->Delete();
		delete treeReader;
	}
	//___________________________________________
	TFile* outputFile = new TFile(("../outputs/" + output + "/" + output + "_MVA_selection.root").c_str(), "recreate");
	outputFile->cd();
	data->Write();
	delete data;
	outputFile->Close();
	delete outputFile;
}

int main(int argc, char *argv[]) { //input, output, N events, truth
	std::map<std::string, std::string> options = getOptions(argc, argv);
	if (options.size() == 0) {
		return 1;
	}
	std::string outputName(options["-o"]);
	makeDirs(outputName);
	//ROOT settings______________________________
	gSystem->Load("libDelphes.so");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetPadGridX(kFALSE);
	gStyle->SetPadGridY(kFALSE);
	//___________________________________________
	if (options["-m"] != "0") printMVASelectionInfo(options["-i"], outputName, options["-t"]); //Output info for MVA-based selection
	//Initialise variables_______________________
	std::cout << "Initialising variables\n";
	int lepton_0, lepton_1, tau_0, tau_1, bJet_0, bJet_1;
	//Low-level variables________________________
	double t_0_pT, t_0_eta, t_0_phi, t_0_mass; //Tau 0 variables
	double t_1_pT, t_1_eta, t_1_phi, t_1_mass; //Tau 1 variables
	double b_0_pT, b_0_eta, b_0_phi, b_0_mass; //b-jet 0 variables
	double b_1_pT, b_1_eta, b_1_phi, b_1_mass; //b-jet 1 variables
	double mPT_pT, mPT_phi; //Missing ET variables
	//___________________________________________
	//Reconstructed variables____________________
	double h_tt_pT, h_tt_eta, h_tt_phi, h_tt_mass; //Higgs 0 variables
	double h_bb_pT, h_bb_eta, h_bb_phi, h_bb_mass; //Higgs 1 variables
	double diH_pT, diH_eta, diH_phi, diH_mass; //di-Higgs variables
	//___________________________________________
	//Global event variables_____________________
	double hT, sT, centrality, eVis; //Global kinematics
	int nJets, nBJets, nTauJets; //Jet multiplicities
	double minJetPT, meanJetPT, maxJetPT; //Global jet pTs
	double minJetMass, meanJetMass, maxJetMass; //Global jet masses
	double minJetEta, meanJetEta, maxJetEta; //Global jet etas
	int nPhotons;
	double sphericityA, spherocityA, aplanarityA, aplanorityA, upsilonA, dShapeA; //Event shapes for all objects
	double sphericityP, spherocityP, aplanarityP, aplanorityP, upsilonP, dShapeP; //Event shapes for primary objects
	//___________________________________________
	//Generator-level variables for regression and cuts
	double gen_t_0_pT, gen_t_0_eta, gen_t_0_phi, gen_t_0_E; //Tau 0 variables
	double gen_t_1_pT, gen_t_1_eta, gen_t_1_phi, gen_t_1_E; //Tau 1 variables
	double gen_b_0_pT, gen_b_0_eta, gen_b_0_phi, gen_b_0_E; //b-jet 0 variables
	double gen_b_1_pT, gen_b_1_eta, gen_b_1_phi, gen_b_1_E; //b-jet 1 variables
	double gen_diH_pT, gen_diH_eta, gen_diH_phi, gen_diH_E, gen_diH_mass; //diHiggs variables
	double gen_h_bb_pT, gen_h_bb_eta, gen_h_bb_phi, gen_h_bb_E; //Higgs->bb variables
	double gen_h_tt_pT, gen_h_tt_eta, gen_h_tt_phi, gen_h_tt_E; //Higgs->tau tau variables
	bool gen_mctMatch; //MC truth match
	//___________________________________________
	double weight; //Event weight
	int nElectrons = 0, nMuons = 0;
	bool eventAccepted = false;
	TTree* e_tau_b_b = new TTree("e_tau_b_b", "e #tau b #bar{b}");
	e_tau_b_b->Branch("t_0_pT", &t_0_pT);
	e_tau_b_b->Branch("t_0_eta", &t_0_eta);
	e_tau_b_b->Branch("t_0_phi", &t_0_phi);
	e_tau_b_b->Branch("t_0_mass", &t_0_mass);
	e_tau_b_b->Branch("t_1_pT", &t_1_pT);
	e_tau_b_b->Branch("t_1_eta", &t_1_eta);
	e_tau_b_b->Branch("t_1_phi", &t_1_phi);
	e_tau_b_b->Branch("t_1_mass", &t_1_mass);
	e_tau_b_b->Branch("b_0_pT", &b_0_pT);
	e_tau_b_b->Branch("b_0_eta", &b_0_eta);
	e_tau_b_b->Branch("b_0_phi", &b_0_phi);
	e_tau_b_b->Branch("b_0_mass", &b_0_mass);
	e_tau_b_b->Branch("b_1_pT", &b_1_pT);
	e_tau_b_b->Branch("b_1_eta", &b_1_eta);
	e_tau_b_b->Branch("b_1_phi", &b_1_phi);
	e_tau_b_b->Branch("b_1_mass", &b_1_mass);
	e_tau_b_b->Branch("mPT_pT", &mPT_pT);
	e_tau_b_b->Branch("mPT_phi", &mPT_phi);
	e_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	e_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	e_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	e_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	e_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	e_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	e_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	e_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	e_tau_b_b->Branch("diH_pT", &diH_pT);
	e_tau_b_b->Branch("diH_eta", &diH_eta);
	e_tau_b_b->Branch("diH_phi", &diH_phi);
	e_tau_b_b->Branch("diH_mass", &diH_mass);
	e_tau_b_b->Branch("hT", &hT);
	e_tau_b_b->Branch("sT", &sT);
	e_tau_b_b->Branch("centrality", &centrality);
	e_tau_b_b->Branch("eVis", &eVis);
	e_tau_b_b->Branch("nJets", &nJets);
	e_tau_b_b->Branch("nBJets", &nBJets);
	e_tau_b_b->Branch("nTauJets", &nTauJets);
	e_tau_b_b->Branch("minJetPT", &minJetPT);
	e_tau_b_b->Branch("meanJetPT", &meanJetPT);
	e_tau_b_b->Branch("maxJetPT", &maxJetPT);
	e_tau_b_b->Branch("minJetMass", &minJetMass);
	e_tau_b_b->Branch("meanJetMass", &meanJetMass);
	e_tau_b_b->Branch("maxJetMass", &maxJetMass);
	e_tau_b_b->Branch("minJetEta", &minJetEta);
	e_tau_b_b->Branch("meanJetEta", &meanJetEta);
	e_tau_b_b->Branch("maxJetEta", &maxJetEta);
	e_tau_b_b->Branch("nPhotons", &nPhotons);
	e_tau_b_b->Branch("sphericityA", &sphericityA);
	e_tau_b_b->Branch("spherocityA", &spherocityA);
	e_tau_b_b->Branch("aplanarityA", &aplanarityA);
	e_tau_b_b->Branch("aplanorityA", &aplanorityA);
	e_tau_b_b->Branch("upsilonA", &upsilonA);
	e_tau_b_b->Branch("dShapeA", &dShapeA);
	e_tau_b_b->Branch("sphericityP", &sphericityP);
	e_tau_b_b->Branch("spherocityP", &spherocityP);
	e_tau_b_b->Branch("aplanarityP", &aplanarityP);
	e_tau_b_b->Branch("aplanorityP", &aplanorityP);
	e_tau_b_b->Branch("upsilonP", &upsilonP);
	e_tau_b_b->Branch("dShapeP", &dShapeP);
	e_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	e_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	e_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	e_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	e_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	e_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	e_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	e_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	e_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	e_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	e_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	e_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	e_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	e_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	e_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	e_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	e_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	e_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	e_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	e_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	e_tau_b_b->Branch("gen_weight", &weight);
	TTree* mu_tau_b_b = new TTree("mu_tau_b_b", "#mu #tau_{h} b #bar{b}");
	mu_tau_b_b->Branch("t_0_pT", &t_0_pT);
	mu_tau_b_b->Branch("t_0_eta", &t_0_eta);
	mu_tau_b_b->Branch("t_0_phi", &t_0_phi);
	mu_tau_b_b->Branch("t_0_mass", &t_0_mass);
	mu_tau_b_b->Branch("t_1_pT", &t_1_pT);
	mu_tau_b_b->Branch("t_1_eta", &t_1_eta);
	mu_tau_b_b->Branch("t_1_phi", &t_1_phi);
	mu_tau_b_b->Branch("t_1_mass", &t_1_mass);
	mu_tau_b_b->Branch("b_0_pT", &b_0_pT);
	mu_tau_b_b->Branch("b_0_eta", &b_0_eta);
	mu_tau_b_b->Branch("b_0_phi", &b_0_phi);
	mu_tau_b_b->Branch("b_0_mass", &b_0_mass);
	mu_tau_b_b->Branch("b_1_pT", &b_1_pT);
	mu_tau_b_b->Branch("b_1_eta", &b_1_eta);
	mu_tau_b_b->Branch("b_1_phi", &b_1_phi);
	mu_tau_b_b->Branch("b_1_mass", &b_1_mass);
	mu_tau_b_b->Branch("mPT_pT", &mPT_pT);
	mu_tau_b_b->Branch("mPT_phi", &mPT_phi);
	mu_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	mu_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	mu_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	mu_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	mu_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	mu_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	mu_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	mu_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	mu_tau_b_b->Branch("diH_pT", &diH_pT);
	mu_tau_b_b->Branch("diH_eta", &diH_eta);
	mu_tau_b_b->Branch("diH_phi", &diH_phi);
	mu_tau_b_b->Branch("diH_mass", &diH_mass);
	mu_tau_b_b->Branch("hT", &hT);
	mu_tau_b_b->Branch("sT", &sT);
	mu_tau_b_b->Branch("centrality", &centrality);
	mu_tau_b_b->Branch("eVis", &eVis);
	mu_tau_b_b->Branch("nJets", &nJets);
	mu_tau_b_b->Branch("nBJets", &nBJets);
	mu_tau_b_b->Branch("nTauJets", &nTauJets);
	mu_tau_b_b->Branch("minJetPT", &minJetPT);
	mu_tau_b_b->Branch("meanJetPT", &meanJetPT);
	mu_tau_b_b->Branch("maxJetPT", &maxJetPT);
	mu_tau_b_b->Branch("minJetMass", &minJetMass);
	mu_tau_b_b->Branch("meanJetMass", &meanJetMass);
	mu_tau_b_b->Branch("maxJetMass", &maxJetMass);
	mu_tau_b_b->Branch("minJetEta", &minJetEta);
	mu_tau_b_b->Branch("meanJetEta", &meanJetEta);
	mu_tau_b_b->Branch("maxJetEta", &maxJetEta);
	mu_tau_b_b->Branch("nPhotons", &nPhotons);
	mu_tau_b_b->Branch("sphericityA", &sphericityA);
	mu_tau_b_b->Branch("spherocityA", &spherocityA);
	mu_tau_b_b->Branch("aplanarityA", &aplanarityA);
	mu_tau_b_b->Branch("aplanorityA", &aplanorityA);
	mu_tau_b_b->Branch("upsilonA", &upsilonA);
	mu_tau_b_b->Branch("dShapeA", &dShapeA);
	mu_tau_b_b->Branch("sphericityP", &sphericityP);
	mu_tau_b_b->Branch("spherocityP", &spherocityP);
	mu_tau_b_b->Branch("aplanarityP", &aplanarityP);
	mu_tau_b_b->Branch("aplanorityP", &aplanorityP);
	mu_tau_b_b->Branch("upsilonP", &upsilonP);
	mu_tau_b_b->Branch("dShapeP", &dShapeP);
	mu_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	mu_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	mu_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	mu_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	mu_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	mu_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	mu_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	mu_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	mu_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	mu_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	mu_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	mu_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	mu_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	mu_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	mu_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	mu_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	mu_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	mu_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	mu_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	mu_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	mu_tau_b_b->Branch("gen_weight", &weight);
	TTree* tau_tau_b_b = new TTree("tau_tau_b_b", "#tau_{h} #tau_{h} b #bar{b}");
	tau_tau_b_b->Branch("t_0_pT", &t_0_pT);
	tau_tau_b_b->Branch("t_0_eta", &t_0_eta);
	tau_tau_b_b->Branch("t_0_phi", &t_0_phi);
	tau_tau_b_b->Branch("t_0_mass", &t_0_mass);
	tau_tau_b_b->Branch("t_1_pT", &t_1_pT);
	tau_tau_b_b->Branch("t_1_eta", &t_1_eta);
	tau_tau_b_b->Branch("t_1_phi", &t_1_phi);
	tau_tau_b_b->Branch("t_1_mass", &t_1_mass);
	tau_tau_b_b->Branch("b_0_pT", &b_0_pT);
	tau_tau_b_b->Branch("b_0_eta", &b_0_eta);
	tau_tau_b_b->Branch("b_0_phi", &b_0_phi);
	tau_tau_b_b->Branch("b_0_mass", &b_0_mass);
	tau_tau_b_b->Branch("b_1_pT", &b_1_pT);
	tau_tau_b_b->Branch("b_1_eta", &b_1_eta);
	tau_tau_b_b->Branch("b_1_phi", &b_1_phi);
	tau_tau_b_b->Branch("b_1_mass", &b_1_mass);
	tau_tau_b_b->Branch("mPT_pT", &mPT_pT);
	tau_tau_b_b->Branch("mPT_phi", &mPT_phi);
	tau_tau_b_b->Branch("h_tt_pT", &h_tt_pT);
	tau_tau_b_b->Branch("h_tt_eta", &h_tt_eta);
	tau_tau_b_b->Branch("h_tt_phi", &h_tt_phi);
	tau_tau_b_b->Branch("h_tt_mass", &h_tt_mass);
	tau_tau_b_b->Branch("h_bb_pT", &h_bb_pT);
	tau_tau_b_b->Branch("h_bb_eta", &h_bb_eta);
	tau_tau_b_b->Branch("h_bb_phi", &h_bb_phi);
	tau_tau_b_b->Branch("h_bb_mass", &h_bb_mass);
	tau_tau_b_b->Branch("diH_pT", &diH_pT);
	tau_tau_b_b->Branch("diH_eta", &diH_eta);
	tau_tau_b_b->Branch("diH_phi", &diH_phi);
	tau_tau_b_b->Branch("diH_mass", &diH_mass);
	tau_tau_b_b->Branch("hT", &hT);
	tau_tau_b_b->Branch("sT", &sT);
	tau_tau_b_b->Branch("centrality", &centrality);
	tau_tau_b_b->Branch("eVis", &eVis);
	tau_tau_b_b->Branch("nJets", &nJets);
	tau_tau_b_b->Branch("nBJets", &nBJets);
	tau_tau_b_b->Branch("nTauJets", &nTauJets);
	tau_tau_b_b->Branch("minJetPT", &minJetPT);
	tau_tau_b_b->Branch("meanJetPT", &meanJetPT);
	tau_tau_b_b->Branch("maxJetPT", &maxJetPT);
	tau_tau_b_b->Branch("minJetMass", &minJetMass);
	tau_tau_b_b->Branch("meanJetMass", &meanJetMass);
	tau_tau_b_b->Branch("maxJetMass", &maxJetMass);
	tau_tau_b_b->Branch("minJetEta", &minJetEta);
	tau_tau_b_b->Branch("meanJetEta", &meanJetEta);
	tau_tau_b_b->Branch("maxJetEta", &maxJetEta);
	tau_tau_b_b->Branch("nPhotons", &nPhotons);
	tau_tau_b_b->Branch("sphericityA", &sphericityA);
	tau_tau_b_b->Branch("spherocityA", &spherocityA);
	tau_tau_b_b->Branch("aplanarityA", &aplanarityA);
	tau_tau_b_b->Branch("aplanorityA", &aplanorityA);
	tau_tau_b_b->Branch("upsilonA", &upsilonA);
	tau_tau_b_b->Branch("dShapeA", &dShapeA);
	tau_tau_b_b->Branch("sphericityP", &sphericityP);
	tau_tau_b_b->Branch("spherocityP", &spherocityP);
	tau_tau_b_b->Branch("aplanarityP", &aplanarityP);
	tau_tau_b_b->Branch("aplanorityP", &aplanorityP);
	tau_tau_b_b->Branch("upsilonP", &upsilonP);
	tau_tau_b_b->Branch("dShapeP", &dShapeP);
	tau_tau_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	tau_tau_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	tau_tau_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	tau_tau_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	tau_tau_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	tau_tau_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	tau_tau_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	tau_tau_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	tau_tau_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	tau_tau_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	tau_tau_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	tau_tau_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	tau_tau_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	tau_tau_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	tau_tau_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	tau_tau_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	tau_tau_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	tau_tau_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	tau_tau_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	tau_tau_b_b->Branch("gen_diH_E", &gen_diH_E);
	tau_tau_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	tau_tau_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	tau_tau_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	tau_tau_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	tau_tau_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	tau_tau_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	tau_tau_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	tau_tau_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	tau_tau_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	tau_tau_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	tau_tau_b_b->Branch("gen_weight", &weight);
	TTree* e_e_b_b = new TTree("e_e_b_b", "e e b #bar{b}");
	e_e_b_b->Branch("t_0_pT", &t_0_pT);
	e_e_b_b->Branch("t_0_eta", &t_0_eta);
	e_e_b_b->Branch("t_0_phi", &t_0_phi);
	e_e_b_b->Branch("t_0_mass", &t_0_mass);
	e_e_b_b->Branch("t_1_pT", &t_1_pT);
	e_e_b_b->Branch("t_1_eta", &t_1_eta);
	e_e_b_b->Branch("t_1_phi", &t_1_phi);
	e_e_b_b->Branch("t_1_mass", &t_1_mass);
	e_e_b_b->Branch("b_0_pT", &b_0_pT);
	e_e_b_b->Branch("b_0_eta", &b_0_eta);
	e_e_b_b->Branch("b_0_phi", &b_0_phi);
	e_e_b_b->Branch("b_0_mass", &b_0_mass);
	e_e_b_b->Branch("b_1_pT", &b_1_pT);
	e_e_b_b->Branch("b_1_eta", &b_1_eta);
	e_e_b_b->Branch("b_1_phi", &b_1_phi);
	e_e_b_b->Branch("b_1_mass", &b_1_mass);
	e_e_b_b->Branch("mPT_pT", &mPT_pT);
	e_e_b_b->Branch("mPT_phi", &mPT_phi);
	e_e_b_b->Branch("h_tt_pT", &h_tt_pT);
	e_e_b_b->Branch("h_tt_eta", &h_tt_eta);
	e_e_b_b->Branch("h_tt_phi", &h_tt_phi);
	e_e_b_b->Branch("h_tt_mass", &h_tt_mass);
	e_e_b_b->Branch("h_bb_pT", &h_bb_pT);
	e_e_b_b->Branch("h_bb_eta", &h_bb_eta);
	e_e_b_b->Branch("h_bb_phi", &h_bb_phi);
	e_e_b_b->Branch("h_bb_mass", &h_bb_mass);
	e_e_b_b->Branch("diH_pT", &diH_pT);
	e_e_b_b->Branch("diH_eta", &diH_eta);
	e_e_b_b->Branch("diH_phi", &diH_phi);
	e_e_b_b->Branch("diH_mass", &diH_mass);
	e_e_b_b->Branch("hT", &hT);
	e_e_b_b->Branch("sT", &sT);
	e_e_b_b->Branch("centrality", &centrality);
	e_e_b_b->Branch("eVis", &eVis);
	e_e_b_b->Branch("nJets", &nJets);
	e_e_b_b->Branch("nBJets", &nBJets);
	e_e_b_b->Branch("nTauJets", &nTauJets);
	e_e_b_b->Branch("minJetPT", &minJetPT);
	e_e_b_b->Branch("meanJetPT", &meanJetPT);
	e_e_b_b->Branch("maxJetPT", &maxJetPT);
	e_e_b_b->Branch("minJetMass", &minJetMass);
	e_e_b_b->Branch("meanJetMass", &meanJetMass);
	e_e_b_b->Branch("maxJetMass", &maxJetMass);
	e_e_b_b->Branch("minJetEta", &minJetEta);
	e_e_b_b->Branch("meanJetEta", &meanJetEta);
	e_e_b_b->Branch("maxJetEta", &maxJetEta);
	e_e_b_b->Branch("nPhotons", &nPhotons);
	e_e_b_b->Branch("sphericityA", &sphericityA);
	e_e_b_b->Branch("spherocityA", &spherocityA);
	e_e_b_b->Branch("aplanarityA", &aplanarityA);
	e_e_b_b->Branch("aplanorityA", &aplanorityA);
	e_e_b_b->Branch("upsilonA", &upsilonA);
	e_e_b_b->Branch("dShapeA", &dShapeA);
	e_e_b_b->Branch("sphericityP", &sphericityP);
	e_e_b_b->Branch("spherocityP", &spherocityP);
	e_e_b_b->Branch("aplanarityP", &aplanarityP);
	e_e_b_b->Branch("aplanorityP", &aplanorityP);
	e_e_b_b->Branch("upsilonP", &upsilonP);
	e_e_b_b->Branch("dShapeP", &dShapeP);
	e_e_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_e_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_e_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_e_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	e_e_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_e_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_e_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_e_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	e_e_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	e_e_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	e_e_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	e_e_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	e_e_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	e_e_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	e_e_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	e_e_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	e_e_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	e_e_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	e_e_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	e_e_b_b->Branch("gen_diH_E", &gen_diH_E);
	e_e_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	e_e_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	e_e_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	e_e_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	e_e_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	e_e_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_e_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_e_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_e_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_e_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	e_e_b_b->Branch("gen_weight", &weight);
	TTree* mu_mu_b_b = new TTree("mu_mu_b_b", "#mu #mu b #bar{b}");
	mu_mu_b_b->Branch("t_0_pT", &t_0_pT);
	mu_mu_b_b->Branch("t_0_eta", &t_0_eta);
	mu_mu_b_b->Branch("t_0_phi", &t_0_phi);
	mu_mu_b_b->Branch("t_0_mass", &t_0_mass);
	mu_mu_b_b->Branch("t_1_pT", &t_1_pT);
	mu_mu_b_b->Branch("t_1_eta", &t_1_eta);
	mu_mu_b_b->Branch("t_1_phi", &t_1_phi);
	mu_mu_b_b->Branch("t_1_mass", &t_1_mass);
	mu_mu_b_b->Branch("b_0_pT", &b_0_pT);
	mu_mu_b_b->Branch("b_0_eta", &b_0_eta);
	mu_mu_b_b->Branch("b_0_phi", &b_0_phi);
	mu_mu_b_b->Branch("b_0_mass", &b_0_mass);
	mu_mu_b_b->Branch("b_1_pT", &b_1_pT);
	mu_mu_b_b->Branch("b_1_eta", &b_1_eta);
	mu_mu_b_b->Branch("b_1_phi", &b_1_phi);
	mu_mu_b_b->Branch("b_1_mass", &b_1_mass);
	mu_mu_b_b->Branch("mPT_pT", &mPT_pT);
	mu_mu_b_b->Branch("mPT_phi", &mPT_phi);
	mu_mu_b_b->Branch("h_tt_pT", &h_tt_pT);
	mu_mu_b_b->Branch("h_tt_eta", &h_tt_eta);
	mu_mu_b_b->Branch("h_tt_phi", &h_tt_phi);
	mu_mu_b_b->Branch("h_tt_mass", &h_tt_mass);
	mu_mu_b_b->Branch("h_bb_pT", &h_bb_pT);
	mu_mu_b_b->Branch("h_bb_eta", &h_bb_eta);
	mu_mu_b_b->Branch("h_bb_phi", &h_bb_phi);
	mu_mu_b_b->Branch("h_bb_mass", &h_bb_mass);
	mu_mu_b_b->Branch("diH_pT", &diH_pT);
	mu_mu_b_b->Branch("diH_eta", &diH_eta);
	mu_mu_b_b->Branch("diH_phi", &diH_phi);
	mu_mu_b_b->Branch("diH_mass", &diH_mass);
	mu_mu_b_b->Branch("hT", &hT);
	mu_mu_b_b->Branch("sT", &sT);
	mu_mu_b_b->Branch("centrality", &centrality);
	mu_mu_b_b->Branch("eVis", &eVis);
	mu_mu_b_b->Branch("nJets", &nJets);
	mu_mu_b_b->Branch("nBJets", &nBJets);
	mu_mu_b_b->Branch("nTauJets", &nTauJets);
	mu_mu_b_b->Branch("minJetPT", &minJetPT);
	mu_mu_b_b->Branch("meanJetPT", &meanJetPT);
	mu_mu_b_b->Branch("maxJetPT", &maxJetPT);
	mu_mu_b_b->Branch("minJetMass", &minJetMass);
	mu_mu_b_b->Branch("meanJetMass", &meanJetMass);
	mu_mu_b_b->Branch("maxJetMass", &maxJetMass);
	mu_mu_b_b->Branch("minJetEta", &minJetEta);
	mu_mu_b_b->Branch("meanJetEta", &meanJetEta);
	mu_mu_b_b->Branch("maxJetEta", &maxJetEta);
	mu_mu_b_b->Branch("nPhotons", &nPhotons);
	mu_mu_b_b->Branch("sphericityA", &sphericityA);
	mu_mu_b_b->Branch("spherocityA", &spherocityA);
	mu_mu_b_b->Branch("aplanarityA", &aplanarityA);
	mu_mu_b_b->Branch("aplanorityA", &aplanorityA);
	mu_mu_b_b->Branch("upsilonA", &upsilonA);
	mu_mu_b_b->Branch("dShapeA", &dShapeA);
	mu_mu_b_b->Branch("sphericityP", &sphericityP);
	mu_mu_b_b->Branch("spherocityP", &spherocityP);
	mu_mu_b_b->Branch("aplanarityP", &aplanarityP);
	mu_mu_b_b->Branch("aplanorityP", &aplanorityP);
	mu_mu_b_b->Branch("upsilonP", &upsilonP);
	mu_mu_b_b->Branch("dShapeP", &dShapeP);
	mu_mu_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_mu_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_mu_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_mu_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	mu_mu_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_mu_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_mu_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_mu_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	mu_mu_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	mu_mu_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	mu_mu_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	mu_mu_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	mu_mu_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	mu_mu_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	mu_mu_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	mu_mu_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	mu_mu_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	mu_mu_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	mu_mu_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	mu_mu_b_b->Branch("gen_diH_E", &gen_diH_E);
	mu_mu_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	mu_mu_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	mu_mu_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	mu_mu_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	mu_mu_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	mu_mu_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_mu_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_mu_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_mu_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_mu_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	mu_mu_b_b->Branch("gen_weight", &weight);
	TTree* e_mu_b_b = new TTree("e_mu_b_b", "e #mu b #bar{b}");
	e_mu_b_b->Branch("t_0_pT", &t_0_pT);
	e_mu_b_b->Branch("t_0_eta", &t_0_eta);
	e_mu_b_b->Branch("t_0_phi", &t_0_phi);
	e_mu_b_b->Branch("t_0_mass", &t_0_mass);
	e_mu_b_b->Branch("t_1_pT", &t_1_pT);
	e_mu_b_b->Branch("t_1_eta", &t_1_eta);
	e_mu_b_b->Branch("t_1_phi", &t_1_phi);
	e_mu_b_b->Branch("t_1_mass", &t_1_mass);
	e_mu_b_b->Branch("b_0_pT", &b_0_pT);
	e_mu_b_b->Branch("b_0_eta", &b_0_eta);
	e_mu_b_b->Branch("b_0_phi", &b_0_phi);
	e_mu_b_b->Branch("b_0_mass", &b_0_mass);
	e_mu_b_b->Branch("b_1_pT", &b_1_pT);
	e_mu_b_b->Branch("b_1_eta", &b_1_eta);
	e_mu_b_b->Branch("b_1_phi", &b_1_phi);
	e_mu_b_b->Branch("b_1_mass", &b_1_mass);
	e_mu_b_b->Branch("mPT_pT", &mPT_pT);
	e_mu_b_b->Branch("mPT_phi", &mPT_phi);
	e_mu_b_b->Branch("h_tt_pT", &h_tt_pT);
	e_mu_b_b->Branch("h_tt_eta", &h_tt_eta);
	e_mu_b_b->Branch("h_tt_phi", &h_tt_phi);
	e_mu_b_b->Branch("h_tt_mass", &h_tt_mass);
	e_mu_b_b->Branch("h_bb_pT", &h_bb_pT);
	e_mu_b_b->Branch("h_bb_eta", &h_bb_eta);
	e_mu_b_b->Branch("h_bb_phi", &h_bb_phi);
	e_mu_b_b->Branch("h_bb_mass", &h_bb_mass);
	e_mu_b_b->Branch("diH_pT", &diH_pT);
	e_mu_b_b->Branch("diH_eta", &diH_eta);
	e_mu_b_b->Branch("diH_phi", &diH_phi);
	e_mu_b_b->Branch("diH_mass", &diH_mass);
	e_mu_b_b->Branch("hT", &hT);
	e_mu_b_b->Branch("sT", &sT);
	e_mu_b_b->Branch("centrality", &centrality);
	e_mu_b_b->Branch("eVis", &eVis);
	e_mu_b_b->Branch("nJets", &nJets);
	e_mu_b_b->Branch("nBJets", &nBJets);
	e_mu_b_b->Branch("nTauJets", &nTauJets);
	e_mu_b_b->Branch("minJetPT", &minJetPT);
	e_mu_b_b->Branch("meanJetPT", &meanJetPT);
	e_mu_b_b->Branch("maxJetPT", &maxJetPT);
	e_mu_b_b->Branch("minJetMass", &minJetMass);
	e_mu_b_b->Branch("meanJetMass", &meanJetMass);
	e_mu_b_b->Branch("maxJetMass", &maxJetMass);
	e_mu_b_b->Branch("minJetEta", &minJetEta);
	e_mu_b_b->Branch("meanJetEta", &meanJetEta);
	e_mu_b_b->Branch("maxJetEta", &maxJetEta);
	e_mu_b_b->Branch("nPhotons", &nPhotons);
	e_mu_b_b->Branch("sphericityA", &sphericityA);
	e_mu_b_b->Branch("spherocityA", &spherocityA);
	e_mu_b_b->Branch("aplanarityA", &aplanarityA);
	e_mu_b_b->Branch("aplanorityA", &aplanorityA);
	e_mu_b_b->Branch("upsilonA", &upsilonA);
	e_mu_b_b->Branch("dShapeA", &dShapeA);
	e_mu_b_b->Branch("sphericityP", &sphericityP);
	e_mu_b_b->Branch("spherocityP", &spherocityP);
	e_mu_b_b->Branch("aplanarityP", &aplanarityP);
	e_mu_b_b->Branch("aplanorityP", &aplanorityP);
	e_mu_b_b->Branch("upsilonP", &upsilonP);
	e_mu_b_b->Branch("dShapeP", &dShapeP);
	e_mu_b_b->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_mu_b_b->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_mu_b_b->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_mu_b_b->Branch("gen_t_0_E", &gen_t_0_E);
	e_mu_b_b->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_mu_b_b->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_mu_b_b->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_mu_b_b->Branch("gen_t_1_E", &gen_t_1_E);
	e_mu_b_b->Branch("gen_b_0_pT", &gen_b_0_pT);
	e_mu_b_b->Branch("gen_b_0_eta", &gen_b_0_eta);
	e_mu_b_b->Branch("gen_b_0_phi", &gen_b_0_phi);
	e_mu_b_b->Branch("gen_b_0_E", &gen_b_0_E);
	e_mu_b_b->Branch("gen_b_1_pT", &gen_b_1_pT);
	e_mu_b_b->Branch("gen_b_1_eta", &gen_b_1_eta);
	e_mu_b_b->Branch("gen_b_1_phi", &gen_b_1_phi);
	e_mu_b_b->Branch("gen_b_1_E", &gen_b_1_E);
	e_mu_b_b->Branch("gen_diH_pT", &gen_diH_pT);
	e_mu_b_b->Branch("gen_diH_eta", &gen_diH_eta);
	e_mu_b_b->Branch("gen_diH_phi", &gen_diH_phi);
	e_mu_b_b->Branch("gen_diH_E", &gen_diH_E);
	e_mu_b_b->Branch("gen_diH_mass", &gen_diH_mass);
	e_mu_b_b->Branch("gen_h_bb_pT", &gen_h_bb_pT);
	e_mu_b_b->Branch("gen_h_bb_eta", &gen_h_bb_eta);
	e_mu_b_b->Branch("gen_h_bb_phi", &gen_h_bb_phi);
	e_mu_b_b->Branch("gen_h_bb_E", &gen_h_bb_E);
	e_mu_b_b->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_mu_b_b->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_mu_b_b->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_mu_b_b->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_mu_b_b->Branch("gen_mctMatch", &gen_mctMatch);
	e_mu_b_b->Branch("gen_weight", &weight);
	std::cout << "Variables initialised\n";
	//___________________________________________
	//Initialise plots___________________________
	std::cout << "Initialising plot\n";
	std::map<std::string, TH1D*> mcTruthPlots;
	TH1D* h_datasetSizes = new TH1D("Dataset_sizes", "Dataset sizes", 7, -0.7, 0.7);
	mcTruthPlots.insert(std::make_pair("cuts", new TH1D("mcTruth_cutFlow", "MC Truth Cuts", 20, -2.0, 2.0)));
	mcTruthPlots.insert(std::make_pair("bMatch", new TH1D("mcTruth_bJetMatching", "#DeltaR(b, jet)", 50, 0.0, 0.5)));
	mcTruthPlots.insert(std::make_pair("tauMatch", new TH1D("mcTruth_tauJetMatching", "#DeltaR(#tau, jet)", 50, 0.0, 0.5)));
	mcTruthPlots.insert(std::make_pair("higgsDecay", new TH1D("mcTruth_higgsDecay", "Higgs product |PID|", 50, 0, 50)));
	TH1D* h_e_tau_b_b_cutFlow;
	TH1D* h_mu_tau_b_b_cutFlow;
	TH1D* h_tau_tau_b_b_cutFlow;
	TH1D* h_e_e_b_b_cutFlow;
	TH1D* h_mu_mu_b_b_cutFlow;
	TH1D* h_e_mu_b_b_cutFlow;
	if (options["-t"] == "0") {
		h_e_tau_b_b_cutFlow = new TH1D("e_tau_b_b_Cut_Flow", "e #tau_{h} b #bar{b} cut flow", 8, -0.8, 0.8);
		h_mu_tau_b_b_cutFlow = new TH1D("mu_tau_b_b_Cut_Flow", "#mu #tau_{h} b #bar{b} cut flow", 8, -0.8, 0.8);
		h_tau_tau_b_b_cutFlow = new TH1D("tau_tau_b_b_Cut_Flow", "#tau_{h} #tau_{h} b #bar{b} cut flow", 7, -0.7, 0.7);
		h_e_e_b_b_cutFlow = new TH1D("e_e_b_b_Cut_Flow", "e e b #bar{b} cut flow", 8, -0.8, 0.8);
		h_mu_mu_b_b_cutFlow = new TH1D("mu_mu_b_b_Cut_Flow", "#mu #mu b #bar{b} cut flow", 8, -0.8, 0.8);
		h_e_mu_b_b_cutFlow = new TH1D("e_mu_b_b_Cut_Flow", "e #mu b #bar{b} cut flow", 8, -0.8, 0.8);
	} else {
		h_e_tau_b_b_cutFlow = new TH1D("e_tau_b_b_Cut_Flow", "e #tau_{h} b #bar{b} cut flow", 9, -0.9, 0.9);
		h_mu_tau_b_b_cutFlow = new TH1D("mu_tau_b_b_Cut_Flow", "#mu #tau_{h} b #bar{b} cut flow", 9, -0.9, 0.9);
		h_tau_tau_b_b_cutFlow = new TH1D("tau_tau_b_b_Cut_Flow", "#tau_{h} #tau_{h} b #bar{b} cut flow", 8, -0.8, 0.8);
		h_e_e_b_b_cutFlow = new TH1D("e_e_b_b_Cut_Flow", "e e b #bar{b} cut flow", 9, -0.9, 0.9);
		h_mu_mu_b_b_cutFlow = new TH1D("mu_mu_b_b_Cut_Flow", "#mu #mu b #bar{b} cut flow", 9, -0.9, 0.9);
		h_e_mu_b_b_cutFlow = new TH1D("e_mu_b_b_Cut_Flow", "e #mu b #bar{b} cut flow", 9, -0.9, 0.9);
	}
	h_datasetSizes->GetXaxis()->SetBinLabel(1, "All");
	h_datasetSizes->GetXaxis()->SetBinLabel(2, "#mu #tau_{h} b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(3, "e #tau_{h} b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(4, "#tau_{h} #tau_{h} b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(5, "e e b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(6, "e #mu b #bar{b}");
	h_datasetSizes->GetXaxis()->SetBinLabel(7, "#mu #mu b #bar{b}");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality b#bar{b}");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "Quality e");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "1 e & 0 #mu");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "OS");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	if (options["-t"] == "1") h_e_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality b#bar{b}");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "Quality #mu");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "1 #mu & 0 e");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "OS");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	if (options["-t"] == "1") h_mu_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau#tau");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "Quality b#bar{b}");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "0 e & 0 #mu");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "OS");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "m_{#tau#tau} Cut");
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{b#bar{b}} Cut");
	if (options["-t"] == "1") h_tau_tau_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "MC truth");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality di-e");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "2 e & 0 #mu");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	if (options["-t"] == "1") h_e_e_b_b_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality e and #mu");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "1 e & 1 #mu");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	if (options["-t"] == "1") h_e_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(2, "Quality di-#mu");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(3, "2 #mu & 0 e");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(5, "Quality b#bar{b}");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(7, "m_{#tau#tau} Cut");
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(8, "m_{b#bar{b}} Cut");
	if (options["-t"] == "1") h_mu_mu_b_b_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(1, "hh->bb#tau#tau check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(2, "hh->bb#tau#tau pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(3, "MC-truth check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(4, "MC-truth pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(5, "b-jets check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(6, "b-jets pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(7, "#taus check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(8, "#taus pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(9, "h->#tau#tau->ee check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(10, "h->#tau#tau->ee pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(11, "h->#tau#tau->e#mu check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(12, "h->#tau#tau->e#mu pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(13, "h->#tau#tau->#mu#mu check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(14, "h->#tau#tau->#mu#mu pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(15, "h->#tau#tau->e#tau_{h} check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(16, "h->#tau#tau->e#tau_{h} pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(17, "h->#tau#tau->#mu#tau_{h} check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(18, "h->#tau#tau->#mu#tau_{h} pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(19, "h->#tau#tau->#tau_{h}#tau_{h} check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(20, "h->#tau#tau->#tau_{h}#tau_{h} pass");
	std::cout << "Plots initialised\n";
	//___________________________________________
	//Load data__________________________________
	if (options["-s"] != "1") return 0; //Run event selection
	std::cout << "Running event selection\n";
	std::vector<std::string> inputs = getInputs(options["-i"]);
	std::cout << inputs.size() << " input files found\n";
	for (int f = 0; f < inputs.size(); f++) {
		std::cout << "Loading file " << f+1 << " of " << inputs.size() << "\n";
		TFile* inputData = TFile::Open(inputs[f].c_str()); //File containing Delphes-simulated MC data
		TTree* eventTree = (TTree*)inputData->Get("Delphes");
		delphesReader* reader = new delphesReader(eventTree);
		std::cout << "Data loaded\n";
		//_______________________________________
		//Loop through events____________________
		Long64_t nEvents = reader->fChain->GetEntriesFast();
		std::cout << "Total number of events in file: " << nEvents << "\n";
		std::vector<int> taus, bJets, electrons, muons;
		TLorentzVector v_tau_0, v_tau_1, v_bJet_0, v_bJet_1, v_higgs_tt, v_higgs_bb, v_diHiggs;
		TLorentzVector v_gen_higgs_bb, v_gen_higgs_tt, v_gen_diHiggs, v_gen_tau_0, v_gen_tau_1, v_gen_bJet_0, v_gen_bJet_1;
		std::cout << "Beginning event loop\n";
		for (Long64_t cEvent = 0; cEvent < nEvents; cEvent++) {
			Long64_t nEvent = reader->LoadTree(cEvent); //Load next event
			reader->fChain->GetEntry(cEvent);
			if (nEvent < 0) break; //Load next event
			if (cEvent % 1000 == 0) std::cout << "Loop: " << cEvent << "/" << nEvents << ", " <<
					100*cEvent/nEvents << "%\n";
			h_datasetSizes->Fill("All", 1);
			eventAccepted = false;
			int hBB = -1, hTauTau = -1;
			/*if (options["-t"] == "1") {
				if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Event is not h->bbtautau
			}*/
			nElectrons = getNElectrons(reader);
			nMuons = getNMuons(reader);
			//Check for mu tau b b finalstates___
			h_mu_tau_b_b_cutFlow->Fill("All", 1);
			electrons.clear();
			muons.clear();
			taus.clear();
			bJets.clear();
			finalstateSet("mu_tau_b_b");
			for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
				if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
						&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
					taus.push_back(i);
				}
				if (reader->Jet_TauTag[i] == 0 && reader->Jet_BTag[i] == 1 && reader->Jet_PT[i] > bJetPTMin
						&& std::abs(reader->Jet_Eta[i]) < bJetEtaMax) { //Quality b jet
					bJets.push_back(i);
				}
			}
			if (taus.size() >= 1) {//Quality tau
				h_mu_tau_b_b_cutFlow->Fill("Quality #tau", 1);
				if (bJets.size() >= 2) {//Quality b jets pairs found
					h_mu_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
					for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
						if (reader->Muon_PT[i] > muPTMin && std::abs(reader->Muon_Eta[i]) < muEtaMax
								&& reader->Muon_IsolationVar[i] < muIsoMax) { //Quality muon
							muons.push_back(i);
						}
					}
					if (muons.size() >= 1) { //Quality muon found
						h_mu_tau_b_b_cutFlow->Fill("Quality #mu", 1);
						if (nMuons == 1 && nElectrons == 0) {
							h_mu_tau_b_b_cutFlow->Fill("1 #mu & 0 e", 1);
							if (getOSTauLeptonPair(reader, &taus, &muons, &tau_0, &lepton_0, false) == true) { //Quality OS pair found
								h_mu_tau_b_b_cutFlow->Fill("OS", 1);
								if (selectBJets(reader, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
									v_tau_1 = getTauLepton(reader, lepton_0, "muon");
									v_tau_0 = getTauHadron(reader, tau_0);
									v_higgs_tt = getHiggs2Taus(reader, v_tau_0, v_tau_1);
									if (!massCut || (v_higgs_tt.M() >= higgsMassMin && v_higgs_tt.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
										h_mu_tau_b_b_cutFlow->Fill("m_{#tau#tau} Cut", 1);
										v_bJet_0 = getBJet(reader, bJet_0);
										v_bJet_1 = getBJet(reader, bJet_1);
										v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
										if (!massCut || (v_higgs_bb.M() >= higgsMassMin && v_higgs_bb.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
											h_mu_tau_b_b_cutFlow->Fill("m_{b#bar{b}} Cut", 1);
											v_diHiggs = getDiHiggs(v_higgs_tt, v_higgs_bb);
											gen_mctMatch = false;
											if (options["-t"] == "1") {
												if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Checks if event is h->bbtautau
												gen_mctMatch = truthCut(inputs[f], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
														tau_0, lepton_0, hBB, hTauTau, "tau:muon",
														&mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
														&v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
											}
											if (options["-t"] == "1" & gen_mctMatch) {
												h_mu_tau_b_b_cutFlow->Fill("MC truth", 1);
											}
											if (debug) std::cout << "Accepted mu_tau_b_b event\n";
											v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
											gen_t_0_pT = v_gen_tau_0.Pt();
											gen_t_0_eta = v_gen_tau_0.Eta();
											gen_t_0_phi = v_gen_tau_0.Phi();
											gen_t_0_E = v_gen_tau_0.E();
											gen_t_1_pT = v_gen_tau_1.Pt();
											gen_t_1_eta = v_gen_tau_1.Eta();
											gen_t_1_phi = v_gen_tau_1.Phi();
											gen_t_1_E = v_gen_tau_1.E();
											gen_b_0_pT = v_gen_bJet_0.Pt();
											gen_b_0_eta = v_gen_bJet_0.Eta();
											gen_b_0_phi = v_gen_bJet_0.Phi();
											gen_b_0_E = v_gen_bJet_0.E();
											gen_b_1_pT = v_gen_bJet_1.Pt();
											gen_b_1_eta = v_gen_bJet_1.Eta();
											gen_b_1_phi = v_gen_bJet_1.Phi();
											gen_b_1_E = v_gen_bJet_1.E();
											gen_diH_pT = v_gen_diHiggs.Pt();
											gen_diH_eta = v_gen_diHiggs.Eta();
											gen_diH_phi = v_gen_diHiggs.Phi();
											gen_diH_E = v_gen_diHiggs.E();
											gen_diH_mass = v_gen_diHiggs.M();
											gen_h_bb_pT = v_gen_higgs_bb.Pt();
											gen_h_bb_eta = v_gen_higgs_bb.Eta();
											gen_h_bb_phi = v_gen_higgs_bb.Phi();
											gen_h_bb_E = v_gen_higgs_bb.E();
											gen_h_tt_pT = v_gen_higgs_tt.Pt();
											gen_h_tt_eta = v_gen_higgs_tt.Eta();
											gen_h_tt_phi = v_gen_higgs_tt.Phi();
											gen_h_tt_E = v_gen_higgs_tt.E();
											t_0_pT = v_tau_0.Pt();
											t_0_eta = v_tau_0.Eta();
											t_0_phi = v_tau_0.Phi();
											t_0_mass = v_tau_0.M();
											t_1_pT = v_tau_1.Pt();
											t_1_eta = v_tau_1.Eta();
											t_1_phi = v_tau_1.Phi();
											t_1_mass = muMass;
											b_0_pT = v_bJet_0.Pt();
											b_0_eta = v_bJet_0.Eta();
											b_0_phi = v_bJet_0.Phi();
											b_0_mass = v_bJet_0.M();
											b_1_pT = v_bJet_1.Pt();
											b_1_eta = v_bJet_1.Eta();
											b_1_phi = v_bJet_1.Phi();
											b_1_mass = v_bJet_1.M();
											mPT_pT = reader->MissingET_MET[0];
											mPT_phi = reader->MissingET_Phi[0];
											h_tt_pT = v_higgs_tt.Pt();
											h_tt_eta = v_higgs_tt.Eta();
											h_tt_phi = v_higgs_tt.Phi();
											h_tt_mass = v_higgs_tt.M();
											h_bb_pT = v_higgs_bb.Pt();
											h_bb_eta = v_higgs_bb.Eta();
											h_bb_phi = v_higgs_bb.Phi();
											h_bb_mass = v_higgs_bb.M();
											diH_pT = v_diHiggs.Pt();
											diH_eta = v_diHiggs.Eta();
											diH_phi = v_diHiggs.Phi();
											diH_mass = v_diHiggs.M();
											getGlobalEventInfo(inputs[f], cEvent,
													&hT, &sT, &centrality, &eVis,
													&nJets, &nBJets, &nTauJets, &nPhotons,
													&minJetPT, &meanJetPT, &maxJetPT,
													&minJetMass, &meanJetMass, &maxJetMass,
													&minJetEta, &meanJetEta, &maxJetEta,
													&sphericityA, &spherocityA,
													&aplanarityA, &aplanorityA,
													&upsilonA, &dShapeA);
											getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
													&sphericityP, &spherocityP,
													&aplanarityP, &aplanorityP,
													&upsilonP, &dShapeP);
											weight = (double)*reader->Event_Weight;
											mu_tau_b_b->Fill();
											h_datasetSizes->Fill("#mu #tau_{h} b #bar{b}", 1);
											eventAccepted = true;
										}
									}
								}
							}
						}
					}
				}
			}
			//___________________________________
			if (eventAccepted) continue;
			//Check for e tau b b finalstates____
			h_e_tau_b_b_cutFlow->Fill("All", 1);
			electrons.clear();
			muons.clear();
			taus.clear();
			bJets.clear();
			finalstateSet("e_tau_b_b");
			for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
				if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
						&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
					taus.push_back(i);
				}
				if (reader->Jet_TauTag[i] == 0 && reader->Jet_BTag[i] == 1 && reader->Jet_PT[i] > bJetPTMin
						&& std::abs(reader->Jet_Eta[i]) < bJetEtaMax) { //Quality b jet
					bJets.push_back(i);
				}
			}
			if (taus.size() >= 1) {//Quality tau
				h_e_tau_b_b_cutFlow->Fill("Quality #tau", 1);
				if (bJets.size() >= 2) {//Quality b jets pairs found
					h_e_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
					for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
						if (reader->Electron_PT[i] > ePTMin && std::abs(reader->Electron_Eta[i]) < eEtaMax
								&& reader->Electron_IsolationVar[i] < eIsoMax) { //Quality electron
							electrons.push_back(i);
						}
					}
					if (electrons.size() >= 1) { //Quality electron found
						h_e_tau_b_b_cutFlow->Fill("Quality e", 1);
						if (nElectrons == 1 && nMuons == 0) {
							h_e_tau_b_b_cutFlow->Fill("1 e & 0 #mu", 1);
							if (getOSTauLeptonPair(reader, &taus, &electrons, &tau_0, &lepton_0, true) == true) { //Quality OS pair found
								h_e_tau_b_b_cutFlow->Fill("OS", 1);
								if (selectBJets(reader, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair selected
									v_tau_1 = getTauLepton(reader, lepton_0, "electron");
									v_tau_0 = getTauHadron(reader, tau_0);
									v_higgs_tt = getHiggs2Taus(reader, v_tau_0, v_tau_1);
									if (!massCut || (v_higgs_tt.M() >= higgsMassMin && v_higgs_tt.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
										h_e_tau_b_b_cutFlow->Fill("m_{#tau#tau} Cut", 1);
										v_bJet_0 = getBJet(reader, bJet_0);
										v_bJet_1 = getBJet(reader, bJet_1);
										v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
										if (!massCut || (v_higgs_bb.M() >= higgsMassMin && v_higgs_bb.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
											h_e_tau_b_b_cutFlow->Fill("m_{b#bar{b}} Cut", 1);
											gen_mctMatch = false;
											if (options["-t"] == "1") {
												if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Checks if event is h->bbtautau
												gen_mctMatch = truthCut(inputs[f], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
														tau_0, lepton_0, hBB, hTauTau, "tau:electron",
														&mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
														&v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
											}
											if (options["-t"] == "1" & gen_mctMatch) {
												h_e_tau_b_b_cutFlow->Fill("MC truth", 1);
											}
											if (debug) std::cout << "Accepted e_tau_b_b event\n";
											v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
											gen_t_0_pT = v_gen_tau_0.Pt();
											gen_t_0_eta = v_gen_tau_0.Eta();
											gen_t_0_phi = v_gen_tau_0.Phi();
											gen_t_0_E = v_gen_tau_0.E();
											gen_t_1_pT = v_gen_tau_1.Pt();
											gen_t_1_eta = v_gen_tau_1.Eta();
											gen_t_1_phi = v_gen_tau_1.Phi();
											gen_t_1_E = v_gen_tau_1.E();
											gen_b_0_pT = v_gen_bJet_0.Pt();
											gen_b_0_eta = v_gen_bJet_0.Eta();
											gen_b_0_phi = v_gen_bJet_0.Phi();
											gen_b_0_E = v_gen_bJet_0.E();
											gen_b_1_pT = v_gen_bJet_1.Pt();
											gen_b_1_eta = v_gen_bJet_1.Eta();
											gen_b_1_phi = v_gen_bJet_1.Phi();
											gen_b_1_E = v_gen_bJet_1.E();
											gen_diH_pT = v_gen_diHiggs.Pt();
											gen_diH_eta = v_gen_diHiggs.Eta();
											gen_diH_phi = v_gen_diHiggs.Phi();
											gen_diH_E = v_gen_diHiggs.E();
											gen_diH_mass = v_gen_diHiggs.M();
											gen_h_bb_pT = v_gen_higgs_bb.Pt();
											gen_h_bb_eta = v_gen_higgs_bb.Eta();
											gen_h_bb_phi = v_gen_higgs_bb.Phi();
											gen_h_bb_E = v_gen_higgs_bb.E();
											gen_h_tt_pT = v_gen_higgs_tt.Pt();
											gen_h_tt_eta = v_gen_higgs_tt.Eta();
											gen_h_tt_phi = v_gen_higgs_tt.Phi();
											gen_h_tt_E = v_gen_higgs_tt.E();
											t_0_pT = v_tau_0.Pt();
											t_0_eta = v_tau_0.Eta();
											t_0_phi = v_tau_0.Phi();
											t_0_mass = v_tau_0.M();
											t_1_pT = v_tau_1.Pt();
											t_1_eta = v_tau_1.Eta();
											t_1_phi = v_tau_1.Phi();
											t_1_mass = eMass;
											b_0_pT = v_bJet_0.Pt();
											b_0_eta = v_bJet_0.Eta();
											b_0_phi = v_bJet_0.Phi();
											b_0_mass = v_bJet_0.M();
											b_1_pT = v_bJet_1.Pt();
											b_1_eta = v_bJet_1.Eta();
											b_1_phi = v_bJet_1.Phi();
											b_1_mass = v_bJet_1.M();
											mPT_pT = reader->MissingET_MET[0];
											mPT_phi = reader->MissingET_Phi[0];
											h_tt_pT = v_higgs_tt.Pt();
											h_tt_eta = v_higgs_tt.Eta();
											h_tt_phi = v_higgs_tt.Phi();
											h_tt_mass = v_higgs_tt.M();
											h_bb_pT = v_higgs_bb.Pt();
											h_bb_eta = v_higgs_bb.Eta();
											h_bb_phi = v_higgs_bb.Phi();
											h_bb_mass = v_higgs_bb.M();
											diH_pT = v_diHiggs.Pt();
											diH_eta = v_diHiggs.Eta();
											diH_phi = v_diHiggs.Phi();
											diH_mass = v_diHiggs.M();
											getGlobalEventInfo(inputs[f], cEvent,
													&hT, &sT, &centrality, &eVis,
													&nJets, &nBJets, &nTauJets, &nPhotons,
													&minJetPT, &meanJetPT, &maxJetPT,
													&minJetMass, &meanJetMass, &maxJetMass,
													&minJetEta, &meanJetEta, &maxJetEta,
													&sphericityA, &spherocityA,
													&aplanarityA, &aplanorityA,
													&upsilonA, &dShapeA);
											getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
													&sphericityP, &spherocityP,
													&aplanarityP, &aplanorityP,
													&upsilonP, &dShapeP);
											weight = (double)*reader->Event_Weight;
											e_tau_b_b->Fill();
											h_datasetSizes->Fill("e #tau_{h} b #bar{b}", 1);
											eventAccepted = true;
										}
									}
								}
							}
						}
					}
				}
			}
			//___________________________________
			if (eventAccepted) continue;
			//Check for tau tau b b finalstates__
			h_tau_tau_b_b_cutFlow->Fill("All", 1);
			electrons.clear();
			muons.clear();
			taus.clear();
			bJets.clear();
			finalstateSet("tau_tau_b_b");
			for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
				if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
						&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
					taus.push_back(i);
				}
				if (reader->Jet_TauTag[i] == 0 && reader->Jet_BTag[i] == 1 && reader->Jet_PT[i] > bJetPTMin
						&& std::abs(reader->Jet_Eta[i]) < bJetEtaMax) { //Quality b jet
					bJets.push_back(i);
				}
			}
			if (taus.size() >= 2) {
				h_tau_tau_b_b_cutFlow->Fill("Quality #tau#tau", 1);
				if (bJets.size() >= 2) {//Quality taus  and b jets pairs found
					h_tau_tau_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
					if (nElectrons == 0 && nMuons == 0) {
						h_tau_tau_b_b_cutFlow->Fill("0 e & 0 #mu", 1);
						if (getOSTauTauPair(reader, &taus, &tau_0, &tau_1) == true) { //Quality OS pair found
							h_tau_tau_b_b_cutFlow->Fill("OS", 1);
							if (selectBJets(reader, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
								v_tau_0 = getTauHadron(reader, tau_0);
								v_tau_1 = getTauHadron(reader, tau_1);
								v_higgs_tt = getHiggs2Taus(reader, v_tau_0, v_tau_1);
								if (!massCut || (v_higgs_tt.M() >= higgsMassMin && v_higgs_tt.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
									h_tau_tau_b_b_cutFlow->Fill("m_{#tau#tau} Cut", 1);
									v_bJet_0 = getBJet(reader, bJet_0);
									v_bJet_1 = getBJet(reader, bJet_1);
									v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
									if (!massCut || (v_higgs_bb.M() >= higgsMassMin && v_higgs_bb.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
										h_tau_tau_b_b_cutFlow->Fill("m_{b#bar{b}} Cut", 1);
										if (options["-t"] == "1") {
											if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Event is not h->bbtautau
										}
										gen_mctMatch = false;
										if (options["-t"] == "1") {
											if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Checks if event is h->bbtautau
											gen_mctMatch = truthCut(inputs[f], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
													tau_0, tau_1, hBB, hTauTau, "tau:tau",
													&mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
													&v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
										}
										if (options["-t"] == "1" & gen_mctMatch) {
											h_tau_tau_b_b_cutFlow->Fill("MC truth", 1);
										}
										if (debug) std::cout << "Accepted tau_tau_b_b event\n";
										v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
										gen_t_0_pT = v_gen_tau_0.Pt();
										gen_t_0_eta = v_gen_tau_0.Eta();
										gen_t_0_phi = v_gen_tau_0.Phi();
										gen_t_0_E = v_gen_tau_0.E();
										gen_t_1_pT = v_gen_tau_1.Pt();
										gen_t_1_eta = v_gen_tau_1.Eta();
										gen_t_1_phi = v_gen_tau_1.Phi();
										gen_t_1_E = v_gen_tau_1.E();
										gen_b_0_pT = v_gen_bJet_0.Pt();
										gen_b_0_eta = v_gen_bJet_0.Eta();
										gen_b_0_phi = v_gen_bJet_0.Phi();
										gen_b_0_E = v_gen_bJet_0.E();
										gen_b_1_pT = v_gen_bJet_1.Pt();
										gen_b_1_eta = v_gen_bJet_1.Eta();
										gen_b_1_phi = v_gen_bJet_1.Phi();
										gen_b_1_E = v_gen_bJet_1.E();
										gen_diH_pT = v_gen_diHiggs.Pt();
										gen_diH_eta = v_gen_diHiggs.Eta();
										gen_diH_phi = v_gen_diHiggs.Phi();
										gen_diH_E = v_gen_diHiggs.E();
										gen_diH_mass = v_gen_diHiggs.M();
										gen_h_bb_pT = v_gen_higgs_bb.Pt();
										gen_h_bb_eta = v_gen_higgs_bb.Eta();
										gen_h_bb_phi = v_gen_higgs_bb.Phi();
										gen_h_bb_E = v_gen_higgs_bb.E();
										gen_h_tt_pT = v_gen_higgs_tt.Pt();
										gen_h_tt_eta = v_gen_higgs_tt.Eta();
										gen_h_tt_phi = v_gen_higgs_tt.Phi();
										gen_h_tt_E = v_gen_higgs_tt.E();
										t_0_pT = v_tau_0.Pt();
										t_0_eta = v_tau_0.Eta();
										t_0_phi = v_tau_0.Phi();
										t_0_mass = v_tau_0.M();
										t_1_pT = v_tau_1.Pt();
										t_1_eta = v_tau_1.Eta();
										t_1_phi = v_tau_1.Phi();
										t_1_mass = v_tau_1.M();
										b_0_pT = v_bJet_0.Pt();
										b_0_eta = v_bJet_0.Eta();
										b_0_phi = v_bJet_0.Phi();
										b_0_mass = v_bJet_0.M();
										b_1_pT = v_bJet_1.Pt();
										b_1_eta = v_bJet_1.Eta();
										b_1_phi = v_bJet_1.Phi();
										b_1_mass = v_bJet_1.M();
										mPT_pT = reader->MissingET_MET[0];
										mPT_phi = reader->MissingET_Phi[0];
										h_tt_pT = v_higgs_tt.Pt();
										h_tt_eta = v_higgs_tt.Eta();
										h_tt_phi = v_higgs_tt.Phi();
										h_tt_mass = v_higgs_tt.M();
										h_bb_pT = v_higgs_bb.Pt();
										h_bb_eta = v_higgs_bb.Eta();
										h_bb_phi = v_higgs_bb.Phi();
										h_bb_mass = v_higgs_bb.M();
										diH_pT = v_diHiggs.Pt();
										diH_eta = v_diHiggs.Eta();
										diH_phi = v_diHiggs.Phi();
										diH_mass = v_diHiggs.M();
										getGlobalEventInfo(inputs[f], cEvent,
												&hT, &sT, &centrality, &eVis,
												&nJets, &nBJets, &nTauJets, &nPhotons,
												&minJetPT, &meanJetPT, &maxJetPT,
												&minJetMass, &meanJetMass, &maxJetMass,
												&minJetEta, &meanJetEta, &maxJetEta,
												&sphericityA, &spherocityA,
												&aplanarityA, &aplanorityA,
												&upsilonA, &dShapeA);
										getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
												&sphericityP, &spherocityP,
												&aplanarityP, &aplanorityP,
												&upsilonP, &dShapeP);
										weight = (double)*reader->Event_Weight;
										tau_tau_b_b->Fill();
										h_datasetSizes->Fill("#tau_{h} #tau_{h} b #bar{b}", 1);
										eventAccepted = true;
									}
								}
							}
						}
					}
				}
			}
			//___________________________________
			if (eventAccepted) continue;
			//Check for mu mu b b finalstates______
			h_mu_mu_b_b_cutFlow->Fill("All", 1);
			electrons.clear();
			muons.clear();
			taus.clear();
			bJets.clear();
			finalstateSet("mu_mu_b_b");
			for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
				if (reader->Muon_PT[i] > muPTMin && std::abs(reader->Muon_Eta[i]) < muEtaMax
						&& reader->Muon_IsolationVar[i] < muIsoMax) { //Quality muons
					muons.push_back(i);
				}
			}
			if (muons.size() >= 2) { //Quality di-muon found
				h_mu_mu_b_b_cutFlow->Fill("Quality di-#mu", 1);
				if (nMuons == 2 && nElectrons == 0) {
					h_mu_mu_b_b_cutFlow->Fill("2 #mu & 0 e", 1);
					if (getOSLeptonLeptonPair(reader, &electrons, &muons, &lepton_0, &lepton_1, "muons") == true) { //Quality OS pair found
						h_mu_mu_b_b_cutFlow->Fill("OS", 1);
						for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
							if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
									&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
								taus.push_back(i);
							}
							if (reader->Jet_TauTag[i] == 0 && reader->Jet_BTag[i] == 1 && reader->Jet_PT[i] > bJetPTMin
									&& std::abs(reader->Jet_Eta[i]) < bJetEtaMax) { //Quality b jet
								bJets.push_back(i);
							}
						}
						if (bJets.size() >= 2) {//Quality b jets pairs found
							if (selectBJets(reader, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
								h_mu_mu_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
								if (taus.size() == 0) { //No quality tau
									h_mu_mu_b_b_cutFlow->Fill("0 #tau_{h}", 1);
									v_tau_0 = getTauLepton(reader, lepton_0, "muon");
									v_tau_1 = getTauLepton(reader, lepton_1, "muon");
									v_higgs_tt = getHiggs2Taus(reader, v_tau_0, v_tau_1);
									if (!massCut || (v_higgs_tt.M() >= higgsMassMin && v_higgs_tt.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
										h_mu_mu_b_b_cutFlow->Fill("m_{#tau#tau} Cut", 1);
										v_bJet_0 = getBJet(reader, bJet_0);
										v_bJet_1 = getBJet(reader, bJet_1);
										v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
										if (!massCut || (v_higgs_bb.M() >= higgsMassMin && v_higgs_bb.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
											h_mu_mu_b_b_cutFlow->Fill("m_{b#bar{b}} Cut", 1);
											gen_mctMatch = false;
											if (options["-t"] == "1") {
												if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Checks if event is h->bbtautau
												gen_mctMatch = truthCut(inputs[f], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
														lepton_0, lepton_1, hBB, hTauTau, "muon:muon",
														&mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
														&v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
											}
											if (options["-t"] == "1" & gen_mctMatch) {
												h_mu_mu_b_b_cutFlow->Fill("MC truth", 1);
											}
											if (debug) std::cout << "Accepted mu_mu_b_b event\n";
											v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
											gen_t_0_pT = v_gen_tau_0.Pt();
											gen_t_0_eta = v_gen_tau_0.Eta();
											gen_t_0_phi = v_gen_tau_0.Phi();
											gen_t_0_E = v_gen_tau_0.E();
											gen_t_1_pT = v_gen_tau_1.Pt();
											gen_t_1_eta = v_gen_tau_1.Eta();
											gen_t_1_phi = v_gen_tau_1.Phi();
											gen_t_1_E = v_gen_tau_1.E();
											gen_b_0_pT = v_gen_bJet_0.Pt();
											gen_b_0_eta = v_gen_bJet_0.Eta();
											gen_b_0_phi = v_gen_bJet_0.Phi();
											gen_b_0_E = v_gen_bJet_0.E();
											gen_b_1_pT = v_gen_bJet_1.Pt();
											gen_b_1_eta = v_gen_bJet_1.Eta();
											gen_b_1_phi = v_gen_bJet_1.Phi();
											gen_b_1_E = v_gen_bJet_1.E();
											gen_diH_pT = v_gen_diHiggs.Pt();
											gen_diH_eta = v_gen_diHiggs.Eta();
											gen_diH_phi = v_gen_diHiggs.Phi();
											gen_diH_E = v_gen_diHiggs.E();
											gen_diH_mass = v_gen_diHiggs.M();
											gen_h_bb_pT = v_gen_higgs_bb.Pt();
											gen_h_bb_eta = v_gen_higgs_bb.Eta();
											gen_h_bb_phi = v_gen_higgs_bb.Phi();
											gen_h_bb_E = v_gen_higgs_bb.E();
											gen_h_tt_pT = v_gen_higgs_tt.Pt();
											gen_h_tt_eta = v_gen_higgs_tt.Eta();
											gen_h_tt_phi = v_gen_higgs_tt.Phi();
											gen_h_tt_E = v_gen_higgs_tt.E();
											t_0_pT = v_tau_0.Pt();
											t_0_eta = v_tau_0.Eta();
											t_0_phi = v_tau_0.Phi();
											t_0_mass = muMass;
											t_1_pT = v_tau_1.Pt();
											t_1_eta = v_tau_1.Eta();
											t_1_phi = v_tau_1.Phi();
											t_1_mass = muMass;
											b_0_pT = v_bJet_0.Pt();
											b_0_eta = v_bJet_0.Eta();
											b_0_phi = v_bJet_0.Phi();
											b_0_mass = v_bJet_0.M();
											b_1_pT = v_bJet_1.Pt();
											b_1_eta = v_bJet_1.Eta();
											b_1_phi = v_bJet_1.Phi();
											b_1_mass = v_bJet_1.M();
											mPT_pT = reader->MissingET_MET[0];
											mPT_phi = reader->MissingET_Phi[0];
											h_tt_pT = v_higgs_tt.Pt();
											h_tt_eta = v_higgs_tt.Eta();
											h_tt_phi = v_higgs_tt.Phi();
											h_tt_mass = v_higgs_tt.M();
											h_bb_pT = v_higgs_bb.Pt();
											h_bb_eta = v_higgs_bb.Eta();
											h_bb_phi = v_higgs_bb.Phi();
											h_bb_mass = v_higgs_bb.M();
											diH_pT = v_diHiggs.Pt();
											diH_eta = v_diHiggs.Eta();
											diH_phi = v_diHiggs.Phi();
											diH_mass = v_diHiggs.M();
											getGlobalEventInfo(inputs[f], cEvent,
													&hT, &sT, &centrality, &eVis,
													&nJets, &nBJets, &nTauJets, &nPhotons,
													&minJetPT, &meanJetPT, &maxJetPT,
													&minJetMass, &meanJetMass, &maxJetMass,
													&minJetEta, &meanJetEta, &maxJetEta,
													&sphericityA, &spherocityA,
													&aplanarityA, &aplanorityA,
													&upsilonA, &dShapeA);
											getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
													&sphericityP, &spherocityP,
													&aplanarityP, &aplanorityP,
													&upsilonP, &dShapeP);
											weight = (double)*reader->Event_Weight;
											mu_mu_b_b->Fill();
											h_datasetSizes->Fill("#mu #mu b #bar{b}", 1);
											eventAccepted = true;
										}
									}
								}
							}
						}
					}
				}
			}
			//___________________________________
			if (eventAccepted) continue;
			//Check for e mu b b finalstates______
			h_e_mu_b_b_cutFlow->Fill("All", 1);
			electrons.clear();
			muons.clear();
			taus.clear();
			bJets.clear();
			finalstateSet("e_mu_b_b");
			for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
				if (reader->Muon_PT[i] > muPTMin && std::abs(reader->Muon_Eta[i]) < muEtaMax
						&& reader->Muon_IsolationVar[i] < muIsoMax) { //Quality muons
					muons.push_back(i);
				}
			}
			for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
				if (reader->Electron_PT[i] > muPTMin && std::abs(reader->Electron_Eta[i]) < muEtaMax
						&& reader->Electron_IsolationVar[i] < muIsoMax) { //Quality electrons
					electrons.push_back(i);
				}
			}
			if (muons.size() >= 1 && electrons.size() >= 1) { //Quality muon and electron found
				h_e_mu_b_b_cutFlow->Fill("Quality e and #mu", 1);
				if (nMuons == 1 && nElectrons == 1) {
					h_e_mu_b_b_cutFlow->Fill("1 e & 1 #mu", 1);
					if (getOSLeptonLeptonPair(reader, &electrons, &muons, &lepton_0, &lepton_1, "mixed") == true) { //Quality OS pair found
						h_e_mu_b_b_cutFlow->Fill("OS", 1);
						for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
							if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
									&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
								taus.push_back(i);
							}
							if (reader->Jet_TauTag[i] == 0 && reader->Jet_BTag[i] == 1 && reader->Jet_PT[i] > bJetPTMin
									&& std::abs(reader->Jet_Eta[i]) < bJetEtaMax) { //Quality b jet
								bJets.push_back(i);
							}
						}
						if (bJets.size() >= 2) {//Quality b jets pairs found
							if (selectBJets(reader, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
								h_e_mu_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
								if (taus.size() == 0) { //No quality tau
									h_e_mu_b_b_cutFlow->Fill("0 #tau_{h}", 1);
									v_tau_0 = getTauLepton(reader, lepton_0, "muon");
									v_tau_1 = getTauLepton(reader, lepton_1, "electron");
									v_higgs_tt = getHiggs2Taus(reader, v_tau_0, v_tau_1);
									if (!massCut || (v_higgs_tt.M() >= higgsMassMin && v_higgs_tt.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
										h_e_mu_b_b_cutFlow->Fill("m_{#tau#tau} Cut", 1);
										v_bJet_0 = getBJet(reader, bJet_0);
										v_bJet_1 = getBJet(reader, bJet_1);
										v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
										if (!massCut || (v_higgs_bb.M() >= higgsMassMin && v_higgs_bb.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
											h_e_mu_b_b_cutFlow->Fill("m_{b#bar{b}} Cut", 1);
											gen_mctMatch = false;
											if (options["-t"] == "1") {
												if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Checks if event is h->bbtautau
												gen_mctMatch = truthCut(inputs[f], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
														lepton_0, lepton_1, hBB, hTauTau, "muon:electron",
														&mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
														&v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
											}
											if (options["-t"] == "1" & gen_mctMatch) {
												h_e_mu_b_b_cutFlow->Fill("MC truth", 1);
											}
											if (debug) std::cout << "Accepted e_mu_b_b event\n";
											v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
											gen_t_0_pT = v_gen_tau_0.Pt();
											gen_t_0_eta = v_gen_tau_0.Eta();
											gen_t_0_phi = v_gen_tau_0.Phi();
											gen_t_0_E = v_gen_tau_0.E();
											gen_t_1_pT = v_gen_tau_1.Pt();
											gen_t_1_eta = v_gen_tau_1.Eta();
											gen_t_1_phi = v_gen_tau_1.Phi();
											gen_t_1_E = v_gen_tau_1.E();
											gen_b_0_pT = v_gen_bJet_0.Pt();
											gen_b_0_eta = v_gen_bJet_0.Eta();
											gen_b_0_phi = v_gen_bJet_0.Phi();
											gen_b_0_E = v_gen_bJet_0.E();
											gen_b_1_pT = v_gen_bJet_1.Pt();
											gen_b_1_eta = v_gen_bJet_1.Eta();
											gen_b_1_phi = v_gen_bJet_1.Phi();
											gen_b_1_E = v_gen_bJet_1.E();
											gen_diH_pT = v_gen_diHiggs.Pt();
											gen_diH_eta = v_gen_diHiggs.Eta();
											gen_diH_phi = v_gen_diHiggs.Phi();
											gen_diH_E = v_gen_diHiggs.E();
											gen_diH_mass = v_gen_diHiggs.M();
											gen_h_bb_pT = v_gen_higgs_bb.Pt();
											gen_h_bb_eta = v_gen_higgs_bb.Eta();
											gen_h_bb_phi = v_gen_higgs_bb.Phi();
											gen_h_bb_E = v_gen_higgs_bb.E();
											gen_h_tt_pT = v_gen_higgs_tt.Pt();
											gen_h_tt_eta = v_gen_higgs_tt.Eta();
											gen_h_tt_phi = v_gen_higgs_tt.Phi();
											gen_h_tt_E = v_gen_higgs_tt.E();
											t_0_pT = v_tau_0.Pt();
											t_0_eta = v_tau_0.Eta();
											t_0_phi = v_tau_0.Phi();
											t_0_mass = muMass;
											t_1_pT = v_tau_1.Pt();
											t_1_eta = v_tau_1.Eta();
											t_1_phi = v_tau_1.Phi();
											t_1_mass = eMass;
											b_0_pT = v_bJet_0.Pt();
											b_0_eta = v_bJet_0.Eta();
											b_0_phi = v_bJet_0.Phi();
											b_0_mass = v_bJet_0.M();
											b_1_pT = v_bJet_1.Pt();
											b_1_eta = v_bJet_1.Eta();
											b_1_phi = v_bJet_1.Phi();
											b_1_mass = v_bJet_1.M();
											mPT_pT = reader->MissingET_MET[0];
											mPT_phi = reader->MissingET_Phi[0];
											h_tt_pT = v_higgs_tt.Pt();
											h_tt_eta = v_higgs_tt.Eta();
											h_tt_phi = v_higgs_tt.Phi();
											h_tt_mass = v_higgs_tt.M();
											h_bb_pT = v_higgs_bb.Pt();
											h_bb_eta = v_higgs_bb.Eta();
											h_bb_phi = v_higgs_bb.Phi();
											h_bb_mass = v_higgs_bb.M();
											diH_pT = v_diHiggs.Pt();
											diH_eta = v_diHiggs.Eta();
											diH_phi = v_diHiggs.Phi();
											diH_mass = v_diHiggs.M();
											getGlobalEventInfo(inputs[f], cEvent,
													&hT, &sT, &centrality, &eVis,
													&nJets, &nBJets, &nTauJets, &nPhotons,
													&minJetPT, &meanJetPT, &maxJetPT,
													&minJetMass, &meanJetMass, &maxJetMass,
													&minJetEta, &meanJetEta, &maxJetEta,
													&sphericityA, &spherocityA,
													&aplanarityA, &aplanorityA,
													&upsilonA, &dShapeA);
											getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
													&sphericityP, &spherocityP,
													&aplanarityP, &aplanorityP,
													&upsilonP, &dShapeP);
											weight = (double)*reader->Event_Weight;
											e_mu_b_b->Fill();
											h_datasetSizes->Fill("e #mu b #bar{b}", 1);
											eventAccepted = true;
										}
									}
								}
							}
						}
					}
				}
			}
			//___________________________________
			if (eventAccepted) continue;
			//Check for e e b b finalstates______
			h_e_e_b_b_cutFlow->Fill("All", 1);
			electrons.clear();
			muons.clear();
			taus.clear();
			bJets.clear();
			finalstateSet("e_W_b_b");
			for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
				if (reader->Electron_PT[i] > muPTMin && std::abs(reader->Electron_Eta[i]) < muEtaMax
						&& reader->Electron_IsolationVar[i] < muIsoMax) { //Quality electrons
					electrons.push_back(i);
				}
			}
			if (electrons.size() >= 2) { //Quality di-electron found
				h_e_e_b_b_cutFlow->Fill("Quality di-e", 1);
				if (nMuons == 0 && nElectrons == 2) {
					h_e_e_b_b_cutFlow->Fill("2 e & 0 #mu", 1);
					if (getOSLeptonLeptonPair(reader, &electrons, &muons, &lepton_0, &lepton_1, "electrons") == true) { //Quality OS pair found
						h_e_e_b_b_cutFlow->Fill("OS", 1);
						for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
							if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
									&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
								taus.push_back(i);
							}
							if (reader->Jet_TauTag[i] == 0 && reader->Jet_BTag[i] == 1 && reader->Jet_PT[i] > bJetPTMin
									&& std::abs(reader->Jet_Eta[i]) < bJetEtaMax) { //Quality b jet
								bJets.push_back(i);
							}
						}
						if (bJets.size() >= 2) {//Quality b jets pairs found
							if (selectBJets(reader, &bJets, &bJet_0, &bJet_1) == true) { //Quality b-jet pair found
								h_e_e_b_b_cutFlow->Fill("Quality b#bar{b}", 1);
								if (taus.size() == 0) { //No quality tau
									h_e_e_b_b_cutFlow->Fill("0 #tau_{h}", 1);
									v_tau_0 = getTauLepton(reader, lepton_0, "electron");
									v_tau_1 = getTauLepton(reader, lepton_1, "electron");
									v_higgs_tt = getHiggs2Taus(reader, v_tau_0, v_tau_1);
									if (!massCut || (v_higgs_tt.M() >= higgsMassMin && v_higgs_tt.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
										h_e_e_b_b_cutFlow->Fill("m_{#tau#tau} Cut", 1);
										v_bJet_0 = getBJet(reader, bJet_0);
										v_bJet_1 = getBJet(reader, bJet_1);
										v_higgs_bb = getHiggs2Bs(v_bJet_0, v_bJet_1);
										if (!massCut || (v_higgs_bb.M() >= higgsMassMin && v_higgs_bb.M() <= higgsMassMax)) { //Reconstructed Higgs pass mass window cut
											h_e_e_b_b_cutFlow->Fill("m_{b#bar{b}} Cut", 1);
											gen_mctMatch = false;
											if (options["-t"] == "1") {
												if (!correctDecayChannel(inputs[f], cEvent, &mcTruthPlots, &hBB, &hTauTau)) continue; //Checks if event is h->bbtautau
												gen_mctMatch = truthCut(inputs[f], cEvent, bJet_0, bJet_1, //Checks final-state selection was correct
														lepton_0, lepton_1, hBB, hTauTau, "electron:electron",
														&mcTruthPlots, &v_gen_higgs_bb, &v_gen_higgs_tt,
														&v_gen_tau_0, &v_gen_tau_1, &v_gen_bJet_0, &v_gen_bJet_1);
											}
											if (options["-t"] == "1" & gen_mctMatch) {
												h_e_e_b_b_cutFlow->Fill("MC truth", 1);
											}
											if (debug) std::cout << "Accepted e_e_b_b event\n";
											v_gen_diHiggs = getDiHiggs(v_gen_higgs_tt, v_gen_higgs_bb);
											gen_t_0_pT = v_gen_tau_0.Pt();
											gen_t_0_eta = v_gen_tau_0.Eta();
											gen_t_0_phi = v_gen_tau_0.Phi();
											gen_t_0_E = v_gen_tau_0.E();
											gen_t_1_pT = v_gen_tau_1.Pt();
											gen_t_1_eta = v_gen_tau_1.Eta();
											gen_t_1_phi = v_gen_tau_1.Phi();
											gen_t_1_E = v_gen_tau_1.E();
											gen_b_0_pT = v_gen_bJet_0.Pt();
											gen_b_0_eta = v_gen_bJet_0.Eta();
											gen_b_0_phi = v_gen_bJet_0.Phi();
											gen_b_0_E = v_gen_bJet_0.E();
											gen_b_1_pT = v_gen_bJet_1.Pt();
											gen_b_1_eta = v_gen_bJet_1.Eta();
											gen_b_1_phi = v_gen_bJet_1.Phi();
											gen_b_1_E = v_gen_bJet_1.E();
											gen_diH_pT = v_gen_diHiggs.Pt();
											gen_diH_eta = v_gen_diHiggs.Eta();
											gen_diH_phi = v_gen_diHiggs.Phi();
											gen_diH_E = v_gen_diHiggs.E();
											gen_diH_mass = v_gen_diHiggs.M();
											gen_h_bb_pT = v_gen_higgs_bb.Pt();
											gen_h_bb_eta = v_gen_higgs_bb.Eta();
											gen_h_bb_phi = v_gen_higgs_bb.Phi();
											gen_h_bb_E = v_gen_higgs_bb.E();
											gen_h_tt_pT = v_gen_higgs_tt.Pt();
											gen_h_tt_eta = v_gen_higgs_tt.Eta();
											gen_h_tt_phi = v_gen_higgs_tt.Phi();
											gen_h_tt_E = v_gen_higgs_tt.E();
											t_0_pT = v_tau_0.Pt();
											t_0_eta = v_tau_0.Eta();
											t_0_phi = v_tau_0.Phi();
											t_0_mass = eMass;
											t_1_pT = v_tau_1.Pt();
											t_1_eta = v_tau_1.Eta();
											t_1_phi = v_tau_1.Phi();
											t_1_mass = eMass;
											b_0_pT = v_bJet_0.Pt();
											b_0_eta = v_bJet_0.Eta();
											b_0_phi = v_bJet_0.Phi();
											b_0_mass = v_bJet_0.M();
											b_1_pT = v_bJet_1.Pt();
											b_1_eta = v_bJet_1.Eta();
											b_1_phi = v_bJet_1.Phi();
											b_1_mass = v_bJet_1.M();
											mPT_pT = reader->MissingET_MET[0];
											mPT_phi = reader->MissingET_Phi[0];
											h_tt_pT = v_higgs_tt.Pt();
											h_tt_eta = v_higgs_tt.Eta();
											h_tt_phi = v_higgs_tt.Phi();
											h_tt_mass = v_higgs_tt.M();
											h_bb_pT = v_higgs_bb.Pt();
											h_bb_eta = v_higgs_bb.Eta();
											h_bb_phi = v_higgs_bb.Phi();
											h_bb_mass = v_higgs_bb.M();
											diH_pT = v_diHiggs.Pt();
											diH_eta = v_diHiggs.Eta();
											diH_phi = v_diHiggs.Phi();
											diH_mass = v_diHiggs.M();
											getGlobalEventInfo(inputs[f], cEvent,
													&hT, &sT, &centrality, &eVis,
													&nJets, &nBJets, &nTauJets, &nPhotons,
													&minJetPT, &meanJetPT, &maxJetPT,
													&minJetMass, &meanJetMass, &maxJetMass,
													&minJetEta, &meanJetEta, &maxJetEta,
													&sphericityA, &spherocityA,
													&aplanarityA, &aplanorityA,
													&upsilonA, &dShapeA);
											getPrimaryEventShapes(v_tau_0, v_tau_1, v_bJet_0, v_bJet_1,
													&sphericityP, &spherocityP,
													&aplanarityP, &aplanorityP,
													&upsilonP, &dShapeP);
											weight = (double)*reader->Event_Weight;
											e_e_b_b->Fill();
											h_datasetSizes->Fill("e e b #bar{b}", 1);
											eventAccepted = true;
										}
									}
								}
							}
						}
					}
				}
			}
			//___________________________________
		}
		std::cout << "Event loop complete\n";
		eventTree->Delete();
		inputData->Close();
		delete reader;
	}
	std::cout << "All files complete\n";
	//___________________________________________
	//Writing plots______________________________
	TFile* outputFile = new TFile(("../outputs/" + outputName + "/" + outputName + ".root").c_str(), "recreate");
	outputFile->cd();
	std::cout << "Creating plots\n";
	TCanvas* c_datasetSizes = new TCanvas();
	c_datasetSizes->SetLogy();
	h_datasetSizes->GetXaxis()->SetTitle("Dataset");
	h_datasetSizes->GetYaxis()->SetTitle("Events");
	h_datasetSizes->Draw();
	h_datasetSizes->Write();
	c_datasetSizes->Print(("../outputs/" + outputName + "/datasetSizes.pdf").c_str());
	delete c_datasetSizes;
	TCanvas* c_e_tau_b_b_cutFlow = new TCanvas();
	c_e_tau_b_b_cutFlow->SetLogy();
	h_e_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_tau_b_b_cutFlow->Draw();
	h_e_tau_b_b_cutFlow->Write();
	c_e_tau_b_b_cutFlow->Print(("../outputs/" + outputName + "/e_tau_b_b_cutFlow.pdf").c_str());
	delete c_e_tau_b_b_cutFlow;
	TCanvas* c_mu_tau_b_b_cutFlow = new TCanvas();
	c_mu_tau_b_b_cutFlow->SetLogy();
	h_mu_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_mu_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_mu_tau_b_b_cutFlow->Draw();
	h_mu_tau_b_b_cutFlow->Write();
	c_mu_tau_b_b_cutFlow->Print(("../outputs/" + outputName + "/mu_tau_b_b_cutFlow.pdf").c_str());
	delete c_mu_tau_b_b_cutFlow;
	TCanvas* c_tau_tau_b_b_cutFlow = new TCanvas();
	c_tau_tau_b_b_cutFlow->SetLogy();
	h_tau_tau_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_tau_tau_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_tau_tau_b_b_cutFlow->Draw();
	h_tau_tau_b_b_cutFlow->Write();
	c_tau_tau_b_b_cutFlow->Print(("../outputs/" + outputName + "/tau_tau_b_b_cutFlow.pdf").c_str());
	delete c_tau_tau_b_b_cutFlow;
	TCanvas* c_mu_mu_b_b_cutFlow = new TCanvas();
	c_mu_mu_b_b_cutFlow->SetLogy();
	h_mu_mu_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_mu_mu_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_mu_mu_b_b_cutFlow->Draw();
	h_mu_mu_b_b_cutFlow->Write();
	c_mu_mu_b_b_cutFlow->Print(("../outputs/" + outputName + "/mu_mu_b_b_cutFlow.pdf").c_str());
	delete c_mu_mu_b_b_cutFlow;
	TCanvas* c_e_mu_b_b_cutFlow = new TCanvas();
	c_e_mu_b_b_cutFlow->SetLogy();
	h_e_mu_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_mu_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_mu_b_b_cutFlow->Draw();
	h_e_mu_b_b_cutFlow->Write();
	c_e_mu_b_b_cutFlow->Print(("../outputs/" + outputName + "/e_mu_b_b_cutFlow.pdf").c_str());
	delete c_e_mu_b_b_cutFlow;
	TCanvas* c_e_e_b_b_cutFlow = new TCanvas();
	c_e_e_b_b_cutFlow->SetLogy();
	h_e_e_b_b_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_e_b_b_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_e_b_b_cutFlow->Draw();
	h_e_e_b_b_cutFlow->Write();
	c_e_e_b_b_cutFlow->Print(("../outputs/" + outputName + "/e_e_b_b_cutFlow.pdf").c_str());
	delete c_e_e_b_b_cutFlow;
	TCanvas* c_mcTruth_cutFlow = new TCanvas();
	mcTruthPlots["cuts"]->GetXaxis()->SetTitle("Cuts");
	mcTruthPlots["cuts"]->GetYaxis()->SetTitle("Events");
	mcTruthPlots["cuts"]->Draw();
	mcTruthPlots["cuts"]->Write();
	c_mcTruth_cutFlow->Print(("../outputs/" + outputName + "/mcTruth_cutFlow.pdf").c_str());
	delete c_mcTruth_cutFlow;
	TCanvas* c_mcTruth_bJetMatch = new TCanvas();
	mcTruthPlots["bMatch"]->GetXaxis()->SetTitle("#DeltaR(b, jet)");
	mcTruthPlots["bMatch"]->GetYaxis()->SetTitle("Events");
	mcTruthPlots["bMatch"]->Draw();
	mcTruthPlots["bMatch"]->Write();
	c_mcTruth_bJetMatch->Print(("../outputs/" + outputName + "/mcTruth_bMatch.pdf").c_str());
	delete c_mcTruth_bJetMatch;
	TCanvas* c_mcTruth_tauJetMatch = new TCanvas();
	mcTruthPlots["tauMatch"]->GetXaxis()->SetTitle("#DeltaR(#tau, jet)");
	mcTruthPlots["tauMatch"]->GetYaxis()->SetTitle("Events");
	mcTruthPlots["tauMatch"]->Draw();
	mcTruthPlots["tauMatch"]->Write();
	c_mcTruth_tauJetMatch->Print(("../outputs/" + outputName + "/mcTruth_tauMatch.pdf").c_str());
	delete c_mcTruth_tauJetMatch;
	TCanvas* c_mcTruth_higgsDecay = new TCanvas();
	mcTruthPlots["higgsDecay"]->GetXaxis()->SetTitle("Higgs product |PID|");
	mcTruthPlots["higgsDecay"]->GetYaxis()->SetTitle("Events");
	mcTruthPlots["higgsDecay"]->Draw();
	mcTruthPlots["higgsDecay"]->Write();
	c_mcTruth_higgsDecay->Print(("../outputs/" + outputName + "/mcTruth_higgsDecay.pdf").c_str());
	delete c_mcTruth_higgsDecay;
	std::cout << "Plots created\n";
	//___________________________________________
	//Save datasets______________________________
	std::cout << "Saving data\n";
	e_tau_b_b->Write();
	mu_tau_b_b->Write();
	tau_tau_b_b->Write();
	mu_mu_b_b->Write();
	e_mu_b_b->Write();
	e_e_b_b->Write();
	std::cout << "Data saved\n";
	outputFile->Close();
	delete outputFile;
	delete e_tau_b_b;
	delete mu_tau_b_b;
	delete tau_tau_b_b;
	delete mu_mu_b_b;
	delete e_mu_b_b;
	delete e_e_b_b;
	delete h_datasetSizes;
	delete h_e_tau_b_b_cutFlow;
	delete h_mu_tau_b_b_cutFlow;
	delete h_tau_tau_b_b_cutFlow;
	delete h_e_e_b_b_cutFlow;
	delete h_mu_mu_b_b_cutFlow;
	delete h_e_mu_b_b_cutFlow;
	//___________________________________________
}
