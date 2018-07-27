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
 	to the selected particles.*/
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

void makeDirs(std::string outputName) {
	/*Makes directory structure for saving outputs*/
	std::vector<std::string> dirs;
	dirs.push_back("../outputs");
	dirs.push_back("../outputs/" + outputName);
	for (std::string dir : dirs) {
		system((char*)("mkdir -p " + dir).c_str());
	}
}

TLorentzVector getHiggs2Taus(TLorentzVector mPT, TLorentzVector t_0, TLorentzVector t_1) {
	/*Returns 4-vector of Higgs->tau tau*/
	TLorentzVector higgs;
	higgs = t_0 + t_1 + mPT;
	return higgs;
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

std::map<std::string, std::string> getOptions(int argc, char* argv[]) {
	/*Interpret input arguments*/
	std::map<std::string, std::string> options;
	options.insert(std::make_pair("-i", "")); //Input mask
	options.insert(std::make_pair("-o", "")); //Output name
	options.insert(std::make_pair("-t", "0")); //Truth info
	options.insert(std::make_pair("-d", "0")); //Debug mode
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

bool getGenHiggs(TClonesArray *branchParticle,
		std::map<std::string, TH1D*>* plots=NULL, int* hTauTau=NULL) {
	/*Point and htautau to the Higgs*/
	if (plots != NULL) (*plots)["cuts"]->Fill("h->#tau#tau check", 1);
	for (int p = 0; p < branchParticle->GetEntriesFast(); ++p) {
		if (std::abs(((GenParticle*)branchParticle->At(p))->PID) == 25) { //Particle is Higgs
			if (((GenParticle*)branchParticle->At(p))->D1 >= 0 && ((GenParticle*)branchParticle->At(p))->D2 >= 0) { //Daughters exists
				if (((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID != 25 &&
						((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID != 25) {
					if (plots != NULL) (*plots)["higgsDecay"]->Fill(std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID));
					if (plots != NULL) (*plots)["higgsDecay"]->Fill(std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID));
					if (std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D1))->PID) == 15
							&& std::abs(((GenParticle*)branchParticle->At(((GenParticle*)branchParticle->At(p))->D2))->PID) == 15) { //Daughters are taus
						if (hTauTau != NULL) *hTauTau = p; //Point to Higgs
						if (plots != NULL) (*plots)["cuts"]->Fill("hh->tau#tau pass", 1);
						return true;
					}
				}
			}
		}
	}
	return false; //h->tautau not found
}

bool truthCut(std::string input, Long64_t cEvent, int l_0, int l_1,
		std::string mode, std::map<std::string, TH1D*>* plots,
		TLorentzVector* v_gen_higgs_tt,
		TLorentzVector* v_gen_tau_0, TLorentzVector* v_gen_tau_1) {
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
	int hTauTau;
	if (!getGenHiggs(branchParticle, plots, &hTauTau)) return false;
	int swap;
	//Check taus_________________________________
	if (debug) std::cout << "Checking taus\n";
	std::vector<std::string> options;
	boost::split(options, mode, boost::is_any_of(":"));
	(*plots)["cuts"]->Fill("#taus check", 1);
	(*plots)["cuts"]->Fill(("h->#tau#tau->" + typeLookup(mode) + " check").c_str(), 1);
	GenParticle *tau_0, *tau_1, *higgs;
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
	*v_gen_higgs_tt = ((GenParticle*)branchParticle->At(hTauTau))->P4();
	*v_gen_tau_0 = tau_0->P4();
	*v_gen_tau_1 = tau_1->P4();
	//___________________________________________
	chain->Delete();
	delete treeReader;
	return true;
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
	//Initialise variables_______________________
	std::cout << "Initialising variables\n";
	int lepton_0, lepton_1, tau_0, tau_1;
	//Low-level variables________________________
	double t_0_pT, t_0_eta, t_0_phi, t_0_mass, t_0_mT; //Tau 0 variables
	double t_1_pT, t_1_eta, t_1_phi, t_1_mass, t_1_mT; //Tau 1 variables
	double mPT_pT, mPT_phi; //Missing ET variables
	//___________________________________________
	//Reconstructed variables____________________
	double h_tt_pT, h_tt_eta, h_tt_phi, h_tt_mass; //Higgs 0 variables
	//___________________________________________
	//Generator-level variables for regression and cuts
	double gen_t_0_pT, gen_t_0_eta, gen_t_0_phi, gen_t_0_E; //Tau 0 variables
	double gen_t_1_pT, gen_t_1_eta, gen_t_1_phi, gen_t_1_E; //Tau 1 variables
	double gen_h_tt_pT, gen_h_tt_eta, gen_h_tt_phi, gen_h_tt_E; //Higgs->tau tau variables
	bool gen_mctMatch; //MC truth match
	//___________________________________________
	double weight; //Event weight
	bool eventAccepted = false;
	TTree* e_tau = new TTree("e_tau", "e #tau");
	e_tau->Branch("t_0_pT", &t_0_pT);
	e_tau->Branch("t_0_eta", &t_0_eta);
	e_tau->Branch("t_0_phi", &t_0_phi);
	e_tau->Branch("t_0_mass", &t_0_mass);
	e_tau->Branch("t_0_mT", &t_0_mT);
	e_tau->Branch("t_1_pT", &t_1_pT);
	e_tau->Branch("t_1_eta", &t_1_eta);
	e_tau->Branch("t_1_phi", &t_1_phi);
	e_tau->Branch("t_1_mass", &t_1_mass);
	e_tau->Branch("t_1_mT", &t_1_mT);
	e_tau->Branch("mPT_pT", &mPT_pT);
	e_tau->Branch("mPT_phi", &mPT_phi);
	e_tau->Branch("h_tt_pT", &h_tt_pT);
	e_tau->Branch("h_tt_eta", &h_tt_eta);
	e_tau->Branch("h_tt_phi", &h_tt_phi);
	e_tau->Branch("h_tt_mass", &h_tt_mass);
	e_tau->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_tau->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_tau->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_tau->Branch("gen_t_0_E", &gen_t_0_E);
	e_tau->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_tau->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_tau->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_tau->Branch("gen_t_1_E", &gen_t_1_E);
	e_tau->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_tau->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_tau->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_tau->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_tau->Branch("gen_mctMatch", &gen_mctMatch);
	e_tau->Branch("gen_weight", &weight);
	TTree* mu_tau = new TTree("mu_tau", "#mu #tau_{h}");
	mu_tau->Branch("t_0_pT", &t_0_pT);
	mu_tau->Branch("t_0_eta", &t_0_eta);
	mu_tau->Branch("t_0_phi", &t_0_phi);
	mu_tau->Branch("t_0_mass", &t_0_mass);
	mu_tau->Branch("t_0_mT", &t_0_mT);
	mu_tau->Branch("t_1_pT", &t_1_pT);
	mu_tau->Branch("t_1_eta", &t_1_eta);
	mu_tau->Branch("t_1_phi", &t_1_phi);
	mu_tau->Branch("t_1_mass", &t_1_mass);
	mu_tau->Branch("t_1_mT", &t_1_mT);
	mu_tau->Branch("mPT_pT", &mPT_pT);
	mu_tau->Branch("mPT_phi", &mPT_phi);
	mu_tau->Branch("h_tt_pT", &h_tt_pT);
	mu_tau->Branch("h_tt_eta", &h_tt_eta);
	mu_tau->Branch("h_tt_phi", &h_tt_phi);
	mu_tau->Branch("h_tt_mass", &h_tt_mass);
	mu_tau->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_tau->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_tau->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_tau->Branch("gen_t_0_E", &gen_t_0_E);
	mu_tau->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_tau->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_tau->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_tau->Branch("gen_t_1_E", &gen_t_1_E);
	mu_tau->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_tau->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_tau->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_tau->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_tau->Branch("gen_mctMatch", &gen_mctMatch);
	mu_tau->Branch("gen_weight", &weight);
	TTree* tau_tau = new TTree("tau_tau", "#tau_{h} #tau_{h}");
	tau_tau->Branch("t_0_pT", &t_0_pT);
	tau_tau->Branch("t_0_eta", &t_0_eta);
	tau_tau->Branch("t_0_phi", &t_0_phi);
	tau_tau->Branch("t_0_mass", &t_0_mass);
	tau_tau->Branch("t_0_mT", &t_0_mT);
	tau_tau->Branch("t_1_pT", &t_1_pT);
	tau_tau->Branch("t_1_eta", &t_1_eta);
	tau_tau->Branch("t_1_phi", &t_1_phi);
	tau_tau->Branch("t_1_mass", &t_1_mass);
	tau_tau->Branch("t_1_mT", &t_1_mT);
	tau_tau->Branch("mPT_pT", &mPT_pT);
	tau_tau->Branch("mPT_phi", &mPT_phi);
	tau_tau->Branch("h_tt_pT", &h_tt_pT);
	tau_tau->Branch("h_tt_eta", &h_tt_eta);
	tau_tau->Branch("h_tt_phi", &h_tt_phi);
	tau_tau->Branch("h_tt_mass", &h_tt_mass);
	tau_tau->Branch("gen_t_0_pT", &gen_t_0_pT);
	tau_tau->Branch("gen_t_0_eta", &gen_t_0_eta);
	tau_tau->Branch("gen_t_0_phi", &gen_t_0_phi);
	tau_tau->Branch("gen_t_0_E", &gen_t_0_E);
	tau_tau->Branch("gen_t_1_pT", &gen_t_1_pT);
	tau_tau->Branch("gen_t_1_eta", &gen_t_1_eta);
	tau_tau->Branch("gen_t_1_phi", &gen_t_1_phi);
	tau_tau->Branch("gen_t_1_E", &gen_t_1_E);
	tau_tau->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	tau_tau->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	tau_tau->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	tau_tau->Branch("gen_h_tt_E", &gen_h_tt_E);
	tau_tau->Branch("gen_mctMatch", &gen_mctMatch);
	tau_tau->Branch("gen_weight", &weight);
	TTree* e_e = new TTree("e_e", "e e");
	e_e->Branch("t_0_pT", &t_0_pT);
	e_e->Branch("t_0_eta", &t_0_eta);
	e_e->Branch("t_0_phi", &t_0_phi);
	e_e->Branch("t_0_mass", &t_0_mass);
	e_e->Branch("t_0_mT", &t_0_mT);
	e_e->Branch("t_1_pT", &t_1_pT);
	e_e->Branch("t_1_eta", &t_1_eta);
	e_e->Branch("t_1_phi", &t_1_phi);
	e_e->Branch("t_1_mass", &t_1_mass);
	e_e->Branch("t_1_mT", &t_1_mT);
	e_e->Branch("mPT_pT", &mPT_pT);
	e_e->Branch("mPT_phi", &mPT_phi);
	e_e->Branch("h_tt_pT", &h_tt_pT);
	e_e->Branch("h_tt_eta", &h_tt_eta);
	e_e->Branch("h_tt_phi", &h_tt_phi);
	e_e->Branch("h_tt_mass", &h_tt_mass);
	e_e->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_e->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_e->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_e->Branch("gen_t_0_E", &gen_t_0_E);
	e_e->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_e->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_e->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_e->Branch("gen_t_1_E", &gen_t_1_E);
	e_e->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_e->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_e->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_e->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_e->Branch("gen_mctMatch", &gen_mctMatch);
	e_e->Branch("gen_weight", &weight);
	TTree* mu_mu = new TTree("mu_mu", "#mu #mu");
	mu_mu->Branch("t_0_pT", &t_0_pT);
	mu_mu->Branch("t_0_eta", &t_0_eta);
	mu_mu->Branch("t_0_phi", &t_0_phi);
	mu_mu->Branch("t_0_mass", &t_0_mass);
	mu_mu->Branch("t_0_mT", &t_0_mT);
	mu_mu->Branch("t_1_pT", &t_1_pT);
	mu_mu->Branch("t_1_eta", &t_1_eta);
	mu_mu->Branch("t_1_phi", &t_1_phi);
	mu_mu->Branch("t_1_mass", &t_1_mass);
	mu_mu->Branch("t_1_mT", &t_1_mT);
	mu_mu->Branch("mPT_pT", &mPT_pT);
	mu_mu->Branch("mPT_phi", &mPT_phi);
	mu_mu->Branch("h_tt_pT", &h_tt_pT);
	mu_mu->Branch("h_tt_eta", &h_tt_eta);
	mu_mu->Branch("h_tt_phi", &h_tt_phi);
	mu_mu->Branch("h_tt_mass", &h_tt_mass);
	mu_mu->Branch("gen_t_0_pT", &gen_t_0_pT);
	mu_mu->Branch("gen_t_0_eta", &gen_t_0_eta);
	mu_mu->Branch("gen_t_0_phi", &gen_t_0_phi);
	mu_mu->Branch("gen_t_0_E", &gen_t_0_E);
	mu_mu->Branch("gen_t_1_pT", &gen_t_1_pT);
	mu_mu->Branch("gen_t_1_eta", &gen_t_1_eta);
	mu_mu->Branch("gen_t_1_phi", &gen_t_1_phi);
	mu_mu->Branch("gen_t_1_E", &gen_t_1_E);
	mu_mu->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	mu_mu->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	mu_mu->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	mu_mu->Branch("gen_h_tt_E", &gen_h_tt_E);
	mu_mu->Branch("gen_mctMatch", &gen_mctMatch);
	mu_mu->Branch("gen_weight", &weight);
	TTree* e_mu = new TTree("e_mu", "e #mu");
	e_mu->Branch("t_0_pT", &t_0_pT);
	e_mu->Branch("t_0_eta", &t_0_eta);
	e_mu->Branch("t_0_phi", &t_0_phi);
	e_mu->Branch("t_0_mass", &t_0_mass);
	e_mu->Branch("t_0_mT", &t_0_mT);
	e_mu->Branch("t_1_pT", &t_1_pT);
	e_mu->Branch("t_1_eta", &t_1_eta);
	e_mu->Branch("t_1_phi", &t_1_phi);
	e_mu->Branch("t_1_mass", &t_1_mass);
	e_mu->Branch("t_1_mT", &t_1_mT);
	e_mu->Branch("mPT_pT", &mPT_pT);
	e_mu->Branch("mPT_phi", &mPT_phi);
	e_mu->Branch("h_tt_pT", &h_tt_pT);
	e_mu->Branch("h_tt_eta", &h_tt_eta);
	e_mu->Branch("h_tt_phi", &h_tt_phi);
	e_mu->Branch("h_tt_mass", &h_tt_mass);
	e_mu->Branch("gen_t_0_pT", &gen_t_0_pT);
	e_mu->Branch("gen_t_0_eta", &gen_t_0_eta);
	e_mu->Branch("gen_t_0_phi", &gen_t_0_phi);
	e_mu->Branch("gen_t_0_E", &gen_t_0_E);
	e_mu->Branch("gen_t_1_pT", &gen_t_1_pT);
	e_mu->Branch("gen_t_1_eta", &gen_t_1_eta);
	e_mu->Branch("gen_t_1_phi", &gen_t_1_phi);
	e_mu->Branch("gen_t_1_E", &gen_t_1_E);
	e_mu->Branch("gen_h_tt_pT", &gen_h_tt_pT);
	e_mu->Branch("gen_h_tt_eta", &gen_h_tt_eta);
	e_mu->Branch("gen_h_tt_phi", &gen_h_tt_phi);
	e_mu->Branch("gen_h_tt_E", &gen_h_tt_E);
	e_mu->Branch("gen_mctMatch", &gen_mctMatch);
	e_mu->Branch("gen_weight", &weight);
	std::cout << "Variables initialised\n";
	//___________________________________________
	//Initialise plots___________________________
	std::cout << "Initialising plot\n";
	std::map<std::string, TH1D*> mcTruthPlots;
	TH1D* h_datasetSizes = new TH1D("Dataset_sizes", "Dataset sizes", 7, -0.7, 0.7);
	mcTruthPlots.insert(std::make_pair("cuts", new TH1D("mcTruth_cutFlow", "MC Truth Cuts", 20, -2.0, 2.0)));
	mcTruthPlots.insert(std::make_pair("tauMatch", new TH1D("mcTruth_tauJetMatching", "#DeltaR(#tau, jet)", 50, 0.0, 0.5)));
	mcTruthPlots.insert(std::make_pair("higgsDecay", new TH1D("mcTruth_higgsDecay", "Higgs product |PID|", 50, 0, 50)));
	TH1D* h_e_tau_cutFlow;
	TH1D* h_mu_tau_cutFlow;
	TH1D* h_tau_tau_cutFlow;
	TH1D* h_e_e_cutFlow;
	TH1D* h_mu_mu_cutFlow;
	TH1D* h_e_mu_cutFlow;
	if (options["-t"] == "0") {
		h_e_tau_cutFlow = new TH1D("e_tau_Cut_Flow", "e #tau_{h}cut flow", 8, -0.8, 0.8);
		h_mu_tau_cutFlow = new TH1D("mu_tau_Cut_Flow", "#mu #tau_{h}cut flow", 8, -0.8, 0.8);
		h_tau_tau_cutFlow = new TH1D("tau_tau_Cut_Flow", "#tau_{h} #tau_{h}cut flow", 7, -0.7, 0.7);
		h_e_e_cutFlow = new TH1D("e_e_Cut_Flow", "e ecut flow", 8, -0.8, 0.8);
		h_mu_mu_cutFlow = new TH1D("mu_mu_Cut_Flow", "#mu #mu cut flow", 8, -0.8, 0.8);
		h_e_mu_cutFlow = new TH1D("e_mu_Cut_Flow", "e #mu cut flow", 8, -0.8, 0.8);
	} else {
		h_e_tau_cutFlow = new TH1D("e_tau_Cut_Flow", "e #tau_{h} cut flow", 9, -0.9, 0.9);
		h_mu_tau_cutFlow = new TH1D("mu_tau_Cut_Flow", "#mu #tau_{h} cut flow", 9, -0.9, 0.9);
		h_tau_tau_cutFlow = new TH1D("tau_tau_Cut_Flow", "#tau_{h} #tau_{h} cut flow", 8, -0.8, 0.8);
		h_e_e_cutFlow = new TH1D("e_e_Cut_Flow", "e e b cut flow", 9, -0.9, 0.9);
		h_mu_mu_cutFlow = new TH1D("mu_mu_Cut_Flow", "#mu #mu cut flow", 9, -0.9, 0.9);
		h_e_mu_cutFlow = new TH1D("e_mu_Cut_Flow", "e #mu cut flow", 9, -0.9, 0.9);
	}
	h_datasetSizes->GetXaxis()->SetBinLabel(1, "All");
	h_datasetSizes->GetXaxis()->SetBinLabel(2, "#mu #tau_{h}");
	h_datasetSizes->GetXaxis()->SetBinLabel(3, "e #tau_{h}");
	h_datasetSizes->GetXaxis()->SetBinLabel(4, "#tau_{h} #tau_{h}");
	h_datasetSizes->GetXaxis()->SetBinLabel(5, "e e");
	h_datasetSizes->GetXaxis()->SetBinLabel(6, "e #mu");
	h_datasetSizes->GetXaxis()->SetBinLabel(7, "#mu #mu");
	h_e_tau_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_tau_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau");
	h_e_tau_cutFlow->GetXaxis()->SetBinLabel(4, "Quality e");
	h_e_tau_cutFlow->GetXaxis()->SetBinLabel(5, "1 e & 0 #mu");
	h_e_tau_cutFlow->GetXaxis()->SetBinLabel(6, "OS");
	if (options["-t"] == "1") h_e_tau_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_mu_tau_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_mu_tau_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau");
	h_mu_tau_cutFlow->GetXaxis()->SetBinLabel(4, "Quality #mu");
	h_mu_tau_cutFlow->GetXaxis()->SetBinLabel(5, "1 #mu & 0 e");
	h_mu_tau_cutFlow->GetXaxis()->SetBinLabel(6, "OS");
	if (options["-t"] == "1") h_mu_tau_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_tau_tau_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_tau_tau_cutFlow->GetXaxis()->SetBinLabel(2, "Quality #tau#tau");
	h_tau_tau_cutFlow->GetXaxis()->SetBinLabel(4, "0 e & 0 #mu");
	h_tau_tau_cutFlow->GetXaxis()->SetBinLabel(5, "OS");
	if (options["-t"] == "1") h_tau_tau_cutFlow->GetXaxis()->SetBinLabel(8, "MC truth");
	h_e_e_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_e_cutFlow->GetXaxis()->SetBinLabel(2, "Quality di-e");
	h_e_e_cutFlow->GetXaxis()->SetBinLabel(3, "2 e & 0 #mu");
	h_e_e_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_e_e_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	if (options["-t"] == "1") h_e_e_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_e_mu_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_e_mu_cutFlow->GetXaxis()->SetBinLabel(2, "Quality e and #mu");
	h_e_mu_cutFlow->GetXaxis()->SetBinLabel(3, "1 e & 1 #mu");
	h_e_mu_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_e_mu_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	if (options["-t"] == "1") h_e_mu_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	h_mu_mu_cutFlow->GetXaxis()->SetBinLabel(1, "All");
	h_mu_mu_cutFlow->GetXaxis()->SetBinLabel(2, "Quality di-#mu");
	h_mu_mu_cutFlow->GetXaxis()->SetBinLabel(3, "2 #mu & 0 e");
	h_mu_mu_cutFlow->GetXaxis()->SetBinLabel(4, "OS");
	h_mu_mu_cutFlow->GetXaxis()->SetBinLabel(6, "0 #tau_{h}");
	if (options["-t"] == "1") h_mu_mu_cutFlow->GetXaxis()->SetBinLabel(9, "MC truth");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(1, "hh->#tau#tau check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(2, "hh->#tau#tau pass");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(3, "MC-truth check");
	mcTruthPlots["cuts"]->GetXaxis()->SetBinLabel(4, "MC-truth pass");
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
	TFile* inputData = TFile::Open(options["-i"].c_str()); //File containing Delphes-simulated MC data
	TTree* eventTree = (TTree*)inputData->Get("Delphes");
	delphesReader* reader = new delphesReader(eventTree);
	std::cout << "Data loaded\n";
	//_______________________________________
	//Loop through events____________________
	Long64_t nEvents = reader->fChain->GetEntriesFast();
	std::cout << "Total number of events in file: " << nEvents << "\n";
	std::vector<int> taus, electrons, muons;
	TLorentzVector v_tau_0, v_tau_1, v_higgs_tt, v_mPT;
	TLorentzVector v_gen_higgs_tt, v_gen_tau_0, v_gen_tau_1;
	bool addMuon, addElectron;
	std::cout << "Beginning event loop\n";
	for (Long64_t cEvent = 0; cEvent < nEvents; cEvent++) {
		if (debug) std::cout << "Loading event " << cEvent << "\n";
		Long64_t nEvent = reader->LoadTree(cEvent); //Load next event
		if (debug) std::cout << "Event loaded, getting data\n";
		reader->fChain->GetEntry(cEvent);
		if (debug) std::cout << "Data loaded\n";
		if (nEvent < 0) break; //Load next event
		if (cEvent % 1000 == 0) std::cout << "Loop: " << cEvent << "/" << nEvents << ", " <<
				100*cEvent/nEvents << "%\n";
		h_datasetSizes->Fill("All", 1);
		eventAccepted = false;
		//Check for mu tau finalstates___
		h_mu_tau_cutFlow->Fill("All", 1);
		electrons.clear();
		muons.clear();
		taus.clear();
		addMuon = false;
		addElectron = false;
		finalstateSet("mu_tau");
		if (debug) std::cout << "Running mu tau\n";
		for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
			if (reader->Muon_PT[i] > muPTMin && std::abs(reader->Muon_Eta[i]) < muEtaMax
					&& reader->Muon_IsolationVar[i] < muIsoMax) { //Quality muon
				muons.push_back(i);
			}
		 	else if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
					&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 1 && !addMuon) { //One quality muon found and no additional muons
			h_mu_tau_cutFlow->Fill("Quality #mu", 1);
			for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
				if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
						&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Additional electron
					addElectron = true;
					break;
				}
			}
			if (!addElectron) { //No additional electrons found
				h_mu_tau_cutFlow->Fill("1 #mu & 0 e", 1);
				for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
					if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
							&& std::abs(reader->Jet_Eta[i]) < tauEtaMax && reader->Jet_Charge[i] != reader->Muon_Charge[muons[0]]) { //Quality OS tau
						taus.push_back(i);
					}
				}
				if (taus.size() >= 1) {//Quality tau
					h_mu_tau_cutFlow->Fill("Quality #tau", 1);
					v_tau_1 = getTauLepton(reader, muons[0], "muon");
					v_tau_0 = getTauHadron(reader, taus[0]);
					v_mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0);
					v_higgs_tt = getHiggs2Taus(v_mPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(options["-i"], cEvent, //Checks final-state selection was correct
								taus[0], muons[0], "tau:muon",
								&mcTruthPlots, &v_gen_higgs_tt,
								&v_gen_tau_0, &v_gen_tau_1);
					}
					if (options["-t"] == "1" & gen_mctMatch) {
						h_mu_tau_cutFlow->Fill("MC truth", 1);
					}
					if (debug) std::cout << "Accepted mu_tau event\n";
					gen_t_0_pT = v_gen_tau_0.Pt();
					gen_t_0_eta = v_gen_tau_0.Eta();
					gen_t_0_phi = v_gen_tau_0.Phi();
					gen_t_0_E = v_gen_tau_0.E();
					gen_t_1_pT = v_gen_tau_1.Pt();
					gen_t_1_eta = v_gen_tau_1.Eta();
					gen_t_1_phi = v_gen_tau_1.Phi();
					gen_t_1_E = v_gen_tau_1.E();
					gen_h_tt_pT = v_gen_higgs_tt.Pt();
					gen_h_tt_eta = v_gen_higgs_tt.Eta();
					gen_h_tt_phi = v_gen_higgs_tt.Phi();
					gen_h_tt_E = v_gen_higgs_tt.E();
					mPT_pT = reader->MissingET_MET[0];
					mPT_phi = reader->MissingET_Phi[0];
					t_0_pT = v_tau_0.Pt();
					t_0_eta = v_tau_0.Eta();
					t_0_phi = v_tau_0.Phi();
					t_0_mass = v_tau_0.M();
					t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, v_mPT));
					t_1_pT = v_tau_1.Pt();
					t_1_eta = v_tau_1.Eta();
					t_1_phi = v_tau_1.Phi();
					t_1_mass = muMass;
					t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, v_mPT));
					h_tt_pT = v_higgs_tt.Pt();
					h_tt_eta = v_higgs_tt.Eta();
					h_tt_phi = v_higgs_tt.Phi();
					h_tt_mass = v_higgs_tt.M();
					weight = (double)*reader->Event_Weight;
					mu_tau->Fill();
					h_datasetSizes->Fill("#mu #tau_{h}", 1);
					eventAccepted = true;
				}
			}
		}
		//___________________________________
		if (eventAccepted) continue;
		//Check for e tau finalstates___
		h_e_tau_cutFlow->Fill("All", 1);
		electrons.clear();
		muons.clear();
		taus.clear();
		addMuon = false;
		addElectron = false;
		finalstateSet("e_tau");
		if (debug) std::cout << "Running e tau\n";
		for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
			if (reader->Electron_PT[i] > ePTMin && std::abs(reader->Electron_Eta[i]) < eEtaMax
					&& reader->Electron_IsolationVar[i] < eIsoMax) { //Quality electron
				electrons.push_back(i);
			}
		 	else if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
					&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Additional electon
				addElectron = true;
				break;
			}
		}
		if (electrons.size() == 1 && !addElectron) { //One quality electon found and no additional electons
			h_e_tau_cutFlow->Fill("Quality e", 1);
			for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
				if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
						&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Additional muons
					addMuon = true;
					break;
				}
			}
			if (!addMuon) { //No additional muons found
				h_e_tau_cutFlow->Fill("1 e & 0 #mu", 1);
				for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
					if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
							&& std::abs(reader->Jet_Eta[i]) < tauEtaMax && reader->Jet_Charge[i] != reader->Electron_Charge[electrons[0]]) { //Quality OS tau
						taus.push_back(i);
					}
				}
				if (taus.size() >= 1) {//Quality tau
					h_e_tau_cutFlow->Fill("Quality #tau", 1);
					v_tau_1 = getTauLepton(reader, electrons[0], "electron");
					v_tau_0 = getTauHadron(reader, taus[0]);
					v_mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0);
					v_higgs_tt = getHiggs2Taus(v_mPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(options["-i"], cEvent, //Checks final-state selection was correct
								taus[0], electrons[0], "tau:electron",
								&mcTruthPlots, &v_gen_higgs_tt,
								&v_gen_tau_0, &v_gen_tau_1);
					}
					if (options["-t"] == "1" & gen_mctMatch) {
						h_e_tau_cutFlow->Fill("MC truth", 1);
					}
					if (debug) std::cout << "Accepted e_tau event\n";
					gen_t_0_pT = v_gen_tau_0.Pt();
					gen_t_0_eta = v_gen_tau_0.Eta();
					gen_t_0_phi = v_gen_tau_0.Phi();
					gen_t_0_E = v_gen_tau_0.E();
					gen_t_1_pT = v_gen_tau_1.Pt();
					gen_t_1_eta = v_gen_tau_1.Eta();
					gen_t_1_phi = v_gen_tau_1.Phi();
					gen_t_1_E = v_gen_tau_1.E();
					gen_h_tt_pT = v_gen_higgs_tt.Pt();
					gen_h_tt_eta = v_gen_higgs_tt.Eta();
					gen_h_tt_phi = v_gen_higgs_tt.Phi();
					gen_h_tt_E = v_gen_higgs_tt.E();
					mPT_pT = reader->MissingET_MET[0];
					mPT_phi = reader->MissingET_Phi[0];
					t_0_pT = v_tau_0.Pt();
					t_0_eta = v_tau_0.Eta();
					t_0_phi = v_tau_0.Phi();
					t_0_mass = v_tau_0.M();
					t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, v_mPT));
					t_1_pT = v_tau_1.Pt();
					t_1_eta = v_tau_1.Eta();
					t_1_phi = v_tau_1.Phi();
					t_1_mass = eMass;
					t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, v_mPT));
					h_tt_pT = v_higgs_tt.Pt();
					h_tt_eta = v_higgs_tt.Eta();
					h_tt_phi = v_higgs_tt.Phi();
					h_tt_mass = v_higgs_tt.M();
					weight = (double)*reader->Event_Weight;
					e_tau->Fill();
					h_datasetSizes->Fill("e #tau_{h}", 1);
					eventAccepted = true;
				}
			}
		}
		//___________________________________
		if (eventAccepted) continue;
		//Check for tau tau finalstates___
		h_tau_tau_cutFlow->Fill("All", 1);
		electrons.clear();
		muons.clear();
		taus.clear();
		addMuon = false;
		addElectron = false;
		finalstateSet("tau_tau");
		if (debug) std::cout << "Running tau tau\n";
		for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
			if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
					&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Additional electon
				addElectron = true;
				break;
			}
		}
		if (!addElectron) { //No additional electons
			for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
				if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
						&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Additional muons
					addMuon = true;
					break;
				}
			}
			if (!addMuon) { //No additional muons found
				h_tau_tau_cutFlow->Fill("0 e & 0 #mu", 1);
				for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
					if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
							&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
						taus.push_back(i);
					}
				}
				if (taus.size() >= 2) {//Quality tau pair
					h_tau_tau_cutFlow->Fill("Quality #tau#tau", 1);
					if (getOSTauTauPair(reader, &taus, &tau_0, &tau_1)) { //OS Tau pair
						h_tau_tau_cutFlow->Fill("OS", 1);
						v_tau_1 = getTauHadron(reader, tau_1);
						v_tau_0 = getTauHadron(reader, tau_0);
						v_mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0);
						v_higgs_tt = getHiggs2Taus(v_mPT, v_tau_0, v_tau_1);
						gen_mctMatch = false;
						if (options["-t"] == "1") {
							gen_mctMatch = truthCut(options["-i"], cEvent, //Checks final-state selection was correct
									tau_0, tau_1, "tau:tau",
									&mcTruthPlots, &v_gen_higgs_tt,
									&v_gen_tau_0, &v_gen_tau_1);
						}
						if (options["-t"] == "1" & gen_mctMatch) {
							h_tau_tau_cutFlow->Fill("MC truth", 1);
						}
						if (debug) std::cout << "Accepted tau_tau event\n";
						gen_t_0_pT = v_gen_tau_0.Pt();
						gen_t_0_eta = v_gen_tau_0.Eta();
						gen_t_0_phi = v_gen_tau_0.Phi();
						gen_t_0_E = v_gen_tau_0.E();
						gen_t_1_pT = v_gen_tau_1.Pt();
						gen_t_1_eta = v_gen_tau_1.Eta();
						gen_t_1_phi = v_gen_tau_1.Phi();
						gen_t_1_E = v_gen_tau_1.E();
						gen_h_tt_pT = v_gen_higgs_tt.Pt();
						gen_h_tt_eta = v_gen_higgs_tt.Eta();
						gen_h_tt_phi = v_gen_higgs_tt.Phi();
						gen_h_tt_E = v_gen_higgs_tt.E();
						mPT_pT = reader->MissingET_MET[0];
						mPT_phi = reader->MissingET_Phi[0];
						t_0_pT = v_tau_0.Pt();
						t_0_eta = v_tau_0.Eta();
						t_0_phi = v_tau_0.Phi();
						t_0_mass = v_tau_0.M();
						t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, v_mPT));
						t_1_pT = v_tau_1.Pt();
						t_1_eta = v_tau_1.Eta();
						t_1_phi = v_tau_1.Phi();
						t_1_mass = v_tau_1.M();
						t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, v_mPT));
						h_tt_pT = v_higgs_tt.Pt();
						h_tt_eta = v_higgs_tt.Eta();
						h_tt_phi = v_higgs_tt.Phi();
						h_tt_mass = v_higgs_tt.M();
						weight = (double)*reader->Event_Weight;
						tau_tau->Fill();
						h_datasetSizes->Fill("#tau_{h} #tau_{h}", 1);
						eventAccepted = true;
					}
				}
			}
		}
		//___________________________________
		if (eventAccepted) continue;
		//Check for mu mu finalstates___
		h_mu_mu_cutFlow->Fill("All", 1);
		electrons.clear();
		muons.clear();
		taus.clear();
		addMuon = false;
		addElectron = false;
		finalstateSet("mu_mu");
		if (debug) std::cout << "Running mu mu\n";
		for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
			if (reader->Muon_PT[i] > muPTMin && std::abs(reader->Muon_Eta[i]) < muEtaMax
					&& reader->Muon_IsolationVar[i] < muIsoMax) { //Quality muon
				muons.push_back(i);
			}
		 	else if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
					&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 2 && !addMuon && 
			reader->Muon_Charge[muons[0]] != reader->Muon_Charge[muons[1]]) { //One quality OS muon pair found and no additional muons
			h_mu_mu_cutFlow->Fill("Quality #mu#mu", 1);
			for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
				if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
						&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Additional electron
					addElectron = true;
					break;
				}
			}
			if (!addElectron) { //No additional electrons found
				h_mu_mu_cutFlow->Fill("2 #mu & 0 e", 1);
				for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
					if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
							&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
						taus.push_back(i);
						break;
					}
				}
				if (taus.size() == 0) {//No taus
					h_mu_mu_cutFlow->Fill("0 #tau", 1);
					v_tau_1 = getTauLepton(reader, muons[1], "muon");
					v_tau_0 = getTauLepton(reader, muons[0], "muon");
					v_mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0);
					v_higgs_tt = getHiggs2Taus(v_mPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(options["-i"], cEvent, //Checks final-state selection was correct
								muons[0], muons[1], "muon:muon",
								&mcTruthPlots, &v_gen_higgs_tt,
								&v_gen_tau_0, &v_gen_tau_1);
					}
					if (options["-t"] == "1" & gen_mctMatch) {
						h_mu_mu_cutFlow->Fill("MC truth", 1);
					}
					if (debug) std::cout << "Accepted mu_mu event\n";
					gen_t_0_pT = v_gen_tau_0.Pt();
					gen_t_0_eta = v_gen_tau_0.Eta();
					gen_t_0_phi = v_gen_tau_0.Phi();
					gen_t_0_E = v_gen_tau_0.E();
					gen_t_1_pT = v_gen_tau_1.Pt();
					gen_t_1_eta = v_gen_tau_1.Eta();
					gen_t_1_phi = v_gen_tau_1.Phi();
					gen_t_1_E = v_gen_tau_1.E();
					gen_h_tt_pT = v_gen_higgs_tt.Pt();
					gen_h_tt_eta = v_gen_higgs_tt.Eta();
					gen_h_tt_phi = v_gen_higgs_tt.Phi();
					gen_h_tt_E = v_gen_higgs_tt.E();
					t_0_pT = v_tau_0.Pt();
					t_0_eta = v_tau_0.Eta();
					t_0_phi = v_tau_0.Phi();
					t_0_mass = muMass;
					mPT_pT = reader->MissingET_MET[0];
					mPT_phi = reader->MissingET_Phi[0];
					t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, v_mPT));
					t_1_pT = v_tau_1.Pt();
					t_1_eta = v_tau_1.Eta();
					t_1_phi = v_tau_1.Phi();
					t_1_mass = muMass;
					t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, v_mPT));
					h_tt_pT = v_higgs_tt.Pt();
					h_tt_eta = v_higgs_tt.Eta();
					h_tt_phi = v_higgs_tt.Phi();
					h_tt_mass = v_higgs_tt.M();
					weight = (double)*reader->Event_Weight;
					mu_mu->Fill();
					h_datasetSizes->Fill("#mu #mu", 1);
					eventAccepted = true;
				}
			}
		}
		//___________________________________
		if (eventAccepted) continue;
		//Check for e mu finalstates___
		h_e_mu_cutFlow->Fill("All", 1);
		electrons.clear();
		muons.clear();
		taus.clear();
		addMuon = false;
		addElectron = false;
		finalstateSet("e_mu");
		if (debug) std::cout << "Running e mu\n";
		for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
			if (reader->Muon_PT[i] > muPTMin && std::abs(reader->Muon_Eta[i]) < muEtaMax
					&& reader->Muon_IsolationVar[i] < muIsoMax) { //Quality muon
				muons.push_back(i);
			}
		 	else if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
					&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 1 && !addMuon) { //One quality muon found and no additional muons
			h_e_mu_cutFlow->Fill("Quality #mu", 1);
			for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
				if (reader->Electron_PT[i] > ePTMin && std::abs(reader->Electron_Eta[i]) < eEtaMax
						&& reader->Electron_IsolationVar[i] < eIsoMax
						&& reader->Electron_Charge[i] != reader->Muon_Charge[muons[0]]) { //Quality OS electron
					electrons.push_back(i);
				}
			 	else if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
						&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Additional electon
					addElectron = true;
					break;
				}
			}
			if (electrons.size() == 1 && !addElectron) { //No additional electrons found
				h_e_mu_cutFlow->Fill("1 #mu & 1 e", 1);
				for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
					if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
							&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
						taus.push_back(i);
						break;
					}
				}
				if (taus.size() == 0) {//No taus
					h_e_mu_cutFlow->Fill("0 #tau", 1);
					v_tau_1 = getTauLepton(reader, electrons[0], "electron");
					v_tau_0 = getTauLepton(reader, muons[0], "muon");
					v_mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0);
					v_higgs_tt = getHiggs2Taus(v_mPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(options["-i"], cEvent, //Checks final-state selection was correct
								muons[0], electrons[0], "muon:electron",
								&mcTruthPlots, &v_gen_higgs_tt,
								&v_gen_tau_0, &v_gen_tau_1);
					}
					if (options["-t"] == "1" & gen_mctMatch) {
						h_e_mu_cutFlow->Fill("MC truth", 1);
					}
					if (debug) std::cout << "Accepted e_mu event\n";
					gen_t_0_pT = v_gen_tau_0.Pt();
					gen_t_0_eta = v_gen_tau_0.Eta();
					gen_t_0_phi = v_gen_tau_0.Phi();
					gen_t_0_E = v_gen_tau_0.E();
					gen_t_1_pT = v_gen_tau_1.Pt();
					gen_t_1_eta = v_gen_tau_1.Eta();
					gen_t_1_phi = v_gen_tau_1.Phi();
					gen_t_1_E = v_gen_tau_1.E();
					gen_h_tt_pT = v_gen_higgs_tt.Pt();
					gen_h_tt_eta = v_gen_higgs_tt.Eta();
					gen_h_tt_phi = v_gen_higgs_tt.Phi();
					gen_h_tt_E = v_gen_higgs_tt.E();
					mPT_pT = reader->MissingET_MET[0];
					mPT_phi = reader->MissingET_Phi[0];
					t_0_pT = v_tau_0.Pt();
					t_0_eta = v_tau_0.Eta();
					t_0_phi = v_tau_0.Phi();
					t_0_mass = muMass;
					t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, v_mPT));
					t_1_pT = v_tau_1.Pt();
					t_1_eta = v_tau_1.Eta();
					t_1_phi = v_tau_1.Phi();
					t_1_mass = muMass;
					t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, v_mPT));
					h_tt_pT = v_higgs_tt.Pt();
					h_tt_eta = v_higgs_tt.Eta();
					h_tt_phi = v_higgs_tt.Phi();
					h_tt_mass = v_higgs_tt.M();
					weight = (double)*reader->Event_Weight;
					mu_mu->Fill();
					h_datasetSizes->Fill("e #mu", 1);
					eventAccepted = true;
				}
			}
		}
		//___________________________________
		if (eventAccepted) continue;
		//Check for e e finalstates___
		h_e_e_cutFlow->Fill("All", 1);
		electrons.clear();
		muons.clear();
		taus.clear();
		addMuon = false;
		addElectron = false;
		finalstateSet("e_e");
		if (debug) std::cout << "Running e e\n";
		for (int i = 0; i < reader->Electron_size; i++) { //Loop through electrons
			if (reader->Electron_PT[i] > ePTMin && std::abs(reader->Electron_Eta[i]) < eEtaMax
					&& reader->Electron_IsolationVar[i] < eIsoMax) { //Quality electron
				electrons.push_back(i);
			}
		 	else if (reader->Electron_PT[i] > ePTMinAdd && std::abs(reader->Electron_Eta[i]) < eEtaMaxAdd
					&& reader->Electron_IsolationVar[i] < eIsoMaxAdd) { //Additional electon
				addElectron = true;
				break;
			}
		}
		if (electrons.size() == 2 && !addElectron && 
			reader->Electron_Charge[electrons[0]] != reader->Electron_Charge[electrons[1]]) { //One quality OS electron pair found and no additional electrons
			h_e_e_cutFlow->Fill("Quality ee", 1);
			for (int i = 0; i < reader->Muon_size; i++) { //Loop through muons
				if (reader->Muon_PT[i] > muPTMinAdd && std::abs(reader->Muon_Eta[i]) < muEtaMaxAdd
						&& reader->Muon_IsolationVar[i] < muIsoMaxAdd) { //Additional muons
					addMuon = true;
					break;
				}
			}
			if (!addMuon) { //No additional muons found
				h_e_e_cutFlow->Fill("2 e & 0 #mu", 1);
				for (int i = 0; i < reader->Jet_size; i++) { //Loop through jets
					if (reader->Jet_TauTag[i] == 1 && reader->Jet_BTag[i] == 0 && reader->Jet_PT[i] > tauPTMin
							&& std::abs(reader->Jet_Eta[i]) < tauEtaMax) { //Quality tau
						taus.push_back(i);
						break;
					}
				}
				if (taus.size() == 0) {//No taus
					h_e_e_cutFlow->Fill("0 #tau", 1);
					v_tau_1 = getTauLepton(reader, electrons[1], "electron");
					v_tau_0 = getTauLepton(reader, electrons[0], "electron");
					v_mPT.SetPtEtaPhiM(reader->MissingET_MET[0], 0.0, reader->MissingET_Phi[0], 0.0);
					v_higgs_tt = getHiggs2Taus(v_mPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(options["-i"], cEvent, //Checks final-state selection was correct
								electrons[0], electrons[1], "electron:electron",
								&mcTruthPlots, &v_gen_higgs_tt,
								&v_gen_tau_0, &v_gen_tau_1);
					}
					if (options["-t"] == "1" & gen_mctMatch) {
						h_e_e_cutFlow->Fill("MC truth", 1);
					}
					if (debug) std::cout << "Accepted mu_mu event\n";
					gen_t_0_pT = v_gen_tau_0.Pt();
					gen_t_0_eta = v_gen_tau_0.Eta();
					gen_t_0_phi = v_gen_tau_0.Phi();
					gen_t_0_E = v_gen_tau_0.E();
					gen_t_1_pT = v_gen_tau_1.Pt();
					gen_t_1_eta = v_gen_tau_1.Eta();
					gen_t_1_phi = v_gen_tau_1.Phi();
					gen_t_1_E = v_gen_tau_1.E();
					gen_h_tt_pT = v_gen_higgs_tt.Pt();
					gen_h_tt_eta = v_gen_higgs_tt.Eta();
					gen_h_tt_phi = v_gen_higgs_tt.Phi();
					gen_h_tt_E = v_gen_higgs_tt.E();
					mPT_pT = reader->MissingET_MET[0];
					mPT_phi = reader->MissingET_Phi[0];
					t_0_pT = v_tau_0.Pt();
					t_0_eta = v_tau_0.Eta();
					t_0_phi = v_tau_0.Phi();
					t_0_mass = eMass;
					t_0_mT = getMT(t_0_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_0, v_mPT));
					t_1_pT = v_tau_1.Pt();
					t_1_eta = v_tau_1.Eta();
					t_1_phi = v_tau_1.Phi();
					t_1_mass = eMass;
					t_1_mT = getMT(t_1_pT, mPT_pT, ROOT::Math::VectorUtil::DeltaPhi(v_tau_1, v_mPT));
					h_tt_pT = v_higgs_tt.Pt();
					h_tt_eta = v_higgs_tt.Eta();
					h_tt_phi = v_higgs_tt.Phi();
					h_tt_mass = v_higgs_tt.M();
					weight = (double)*reader->Event_Weight;
					e_e->Fill();
					h_datasetSizes->Fill("e e", 1);
					eventAccepted = true;
				}
			}
		}
		//___________________________________
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
	TCanvas* c_e_tau_cutFlow = new TCanvas();
	c_e_tau_cutFlow->SetLogy();
	h_e_tau_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_tau_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_tau_cutFlow->Draw();
	h_e_tau_cutFlow->Write();
	c_e_tau_cutFlow->Print(("../outputs/" + outputName + "/e_tau_cutFlow.pdf").c_str());
	delete c_e_tau_cutFlow;
	TCanvas* c_mu_tau_cutFlow = new TCanvas();
	c_mu_tau_cutFlow->SetLogy();
	h_mu_tau_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_mu_tau_cutFlow->GetYaxis()->SetTitle("Events");
	h_mu_tau_cutFlow->Draw();
	h_mu_tau_cutFlow->Write();
	c_mu_tau_cutFlow->Print(("../outputs/" + outputName + "/mu_tau_cutFlow.pdf").c_str());
	delete c_mu_tau_cutFlow;
	TCanvas* c_tau_tau_cutFlow = new TCanvas();
	c_tau_tau_cutFlow->SetLogy();
	h_tau_tau_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_tau_tau_cutFlow->GetYaxis()->SetTitle("Events");
	h_tau_tau_cutFlow->Draw();
	h_tau_tau_cutFlow->Write();
	c_tau_tau_cutFlow->Print(("../outputs/" + outputName + "/tau_tau_cutFlow.pdf").c_str());
	delete c_tau_tau_cutFlow;
	TCanvas* c_mu_mu_cutFlow = new TCanvas();
	c_mu_mu_cutFlow->SetLogy();
	h_mu_mu_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_mu_mu_cutFlow->GetYaxis()->SetTitle("Events");
	h_mu_mu_cutFlow->Draw();
	h_mu_mu_cutFlow->Write();
	c_mu_mu_cutFlow->Print(("../outputs/" + outputName + "/mu_mu_cutFlow.pdf").c_str());
	delete c_mu_mu_cutFlow;
	TCanvas* c_e_mu_cutFlow = new TCanvas();
	c_e_mu_cutFlow->SetLogy();
	h_e_mu_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_mu_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_mu_cutFlow->Draw();
	h_e_mu_cutFlow->Write();
	c_e_mu_cutFlow->Print(("../outputs/" + outputName + "/e_mu_cutFlow.pdf").c_str());
	delete c_e_mu_cutFlow;
	TCanvas* c_e_e_cutFlow = new TCanvas();
	c_e_e_cutFlow->SetLogy();
	h_e_e_cutFlow->GetXaxis()->SetTitle("Cuts");
	h_e_e_cutFlow->GetYaxis()->SetTitle("Events");
	h_e_e_cutFlow->Draw();
	h_e_e_cutFlow->Write();
	c_e_e_cutFlow->Print(("../outputs/" + outputName + "/e_e_cutFlow.pdf").c_str());
	delete c_e_e_cutFlow;
	TCanvas* c_mcTruth_cutFlow = new TCanvas();
	mcTruthPlots["cuts"]->GetXaxis()->SetTitle("Cuts");
	mcTruthPlots["cuts"]->GetYaxis()->SetTitle("Events");
	mcTruthPlots["cuts"]->Draw();
	mcTruthPlots["cuts"]->Write();
	c_mcTruth_cutFlow->Print(("../outputs/" + outputName + "/mcTruth_cutFlow.pdf").c_str());
	delete c_mcTruth_cutFlow;
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
	e_tau->Write();
	mu_tau->Write();
	tau_tau->Write();
	mu_mu->Write();
	e_mu->Write();
	e_e->Write();
	std::cout << "Data saved\n";
	outputFile->Close();
	delete outputFile;
	delete e_tau;
	delete mu_tau;
	delete tau_tau;
	delete mu_mu;
	delete e_mu;
	delete e_e;
	delete h_datasetSizes;
	delete h_e_tau_cutFlow;
	delete h_mu_tau_cutFlow;
	delete h_tau_tau_cutFlow;
	delete h_e_e_cutFlow;
	delete h_mu_mu_cutFlow;
	delete h_e_mu_cutFlow;
	//___________________________________________
}