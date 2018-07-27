/*
 * main.cc
 *
 *  Created on: 8 Apr 2016
 *      Author: giles giles.strong@outlook.com
 */
//Local
#include "main.hh"

bool debug = false;

bool getOSTauTauPair(TClonesArray* jets, std::vector<int>* taus, int* tau_0, int* tau_1) {
	/*Checks whether an OS tau-tau pair exists, and if so returns true and points tau and lepton to
	the	selected particles*/
	Jet *tau0, *tau1;
	std::vector<std::pair<int, int> > pairs; //Initialise array for OS tau pairs
	if (debug) std::cout << taus->size() << " tau jets found\n";
	for (int t0 : *taus) { //Loop through taus
		for (int t1 : *taus) {
			if (t0 == t1) continue;
			tau0 = (Jet*)jets->At(t0);
			tau1 = (Jet*)jets->At(t1);
			if (tau0->Charge != tau1->Charge) {
				pairs.push_back(std::make_pair(t0, t1));
			}
		}
	}
	if (debug) std::cout << pairs.size() << " OS tau-jet pairs found\n";
	if (pairs.size() == 1) { //Only one OS pair found
		tau0 = (Jet*)jets->At(pairs[0].first);
		tau1 = (Jet*)jets->At(pairs[0].second);
		if (tau0->PT < tau1->PT) { //Order taus by pT
			*tau_1 = pairs[0].first;
			*tau_0 = pairs[0].second;
		}
		return true;
	} else if (pairs.size() > 1) { //Multiple OS pairs: select highest summed |pT| pair
		double pT, highestPT = 0;
		std::pair<int, int> best;
		for (std::pair<int, int> p : pairs) { //Loop through pairs
			tau0 = (Jet*)jets->At(p.first);
			tau1 = (Jet*)jets->At(p.second);
			pT = tau0->PT+tau1->PT;
			if (pT > highestPT) {
				best = p;
				highestPT = pT;
			}
		}
		*tau_0 = best.first;
		*tau_1 = best.second;
		tau0 = (Jet*)jets->At(best.first);
		tau1 = (Jet*)jets->At(best.second);
		if (tau0->PT < tau1->PT) { //Order taus by pT
			*tau_1 = best.first;
			*tau_0 = best.second;
		}
		return true;
	} else { //No OS pairs found
		return false;
	}
	return false;
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

TLorentzVector getHiggs2Taus(MissingET* mpt, TLorentzVector t_0, TLorentzVector t_1) {
	/*Returns 4-vector of Higgs->tau tau*/
	TLorentzVector higgs, mPT;
	mPT.SetPtEtaPhiM(mpt->MET, 0.0, mpt->Phi, 0.0); //TODO Check this
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

bool truthCut(TClonesArray *branchParticle, TClonesArray *branchElectron,
			  TClonesArray *branchMuon, TClonesArray *branchJet,
			  int l_0, int l_1,
			  std::string mode, std::map<std::string, TH1D*>* plots,
		      TLorentzVector* v_gen_higgs_tt,
			  TLorentzVector* v_gen_tau_0, TLorentzVector* v_gen_tau_1) {
	/*Checks whether selected final states are correct*/
	double jetRadius = 0.5;
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
	std::cout << "Running event selection\n";
	TChain *chain = new TChain("Delphes");
	chain->Add(options["-i"].c_str());
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("MuonLoose");
	TClonesArray *branchJet = treeReader->UseBranch("JetPUPPI");
	TClonesArray *branchMissingET = treeReader->UseBranch("PuppiMissingET");
	TClonesArray *branchWeights = treeReader->UseBranch("Weight");
	std::cout << "Data loaded\n";
	//_______________________________________
	//Loop through events____________________
	Long64_t nEvents = treeReader->GetEntries();
	std::cout << "Total number of events in file: " << nEvents << "\n";
	std::vector<int> taus, bJets, electrons, muons;
	bool addMuon, addElectron;
	TLorentzVector v_tau_0, v_tau_1, v_bJet_0;
	Jet* tmpJet;
	Electron* tmpElectron;
	Muon* tmpMuon;
	MissingET* tmpMPT;
	Weight* tmpWeight;
	std::cout << "Beginning event loop\n";
	for (Long64_t cEvent = 0; cEvent < nEvents; cEvent++) {
		if (debug) std::cout << "Loading event " << cEvent << "\n";
		treeReader->ReadEntry(cEvent); //Load next event
		if (debug) std::cout << "Event loaded, getting data\n";
		if (debug) std::cout << "Data loaded\n";
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
		for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
			tmpMuon = (Muon*) branchMuon->At(i);
			if (tmpMuon->PT > muPTMin && std::abs(tmpMuon->Eta) < muEtaMax
					&& tmpMuon->IsolationVar < muIsoMax) { //Quality muon
				muons.push_back(i);
			} else if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
					&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 1 && !addMuon) { //One quality muon found and no additional muons
			h_mu_tau_cutFlow->Fill("Quality #mu", 1);
			tmpMuon = (Muon*) branchMuon->At(muons[0]);
			for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through electrons
				tmpElectron = (Electron*) branchElectron->At(i);
				if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
						&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electron
					addElectron = true;
					break;
				}
			}
			if (!addElectron) { //No additional electrons found
				h_mu_tau_cutFlow->Fill("1 #mu & 0 e", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax
							&& tmpJet->Charge != tmpMuon->Charge) { //Quality  OS tau
						taus.push_back(i);
					}
				}
				if (taus.size() >= 1) {//Quality tau
					h_mu_tau_cutFlow->Fill("Quality #tau", 1);
					v_tau_1 = tmpMuon->P4()
					tmpJet = (Jet*)branchJet->At(taus[0]);
					v_tau_0 = tmpJet->P4();
					tmpMPT = (MissingET*)branchMissingET->At(0);
					v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(branchParticle, branchElectron,
								branchMuon, branchJet,
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
					mPT_pT = tmpMPT->MET;
					mPT_phi = tmpMPT->Phi;
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
					tmpWeight = (Weight*)branchWeights->At(0);
					weight = tmpWeight->Weight;
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
		for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through Electons
			tmpElectron = (Electron*) branchElectron->At(i);
			if (tmpElectron->PT > ePTMin && std::abs(tmpElectron->Eta) < eEtaMax
					&& tmpElectron->IsolationVar < eIsoMax) { //Quality Electron
				electrons.push_back(i);
			} else if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
					&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electon
				addElectron = true;
				break;
			}
		}
		if (electrons.size() == 1 && !addElectron) { //One quality electon found and no additional electons
			h_e_tau_cutFlow->Fill("Quality e", 1);
			tmpElectron = (Electron*) branchElectron->At(electrons[0]);
			for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
				tmpMuon = (Muon*) branchMuon->At(i);
				if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
						&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
					addMuon = true;
					break;
				}
			}
			if (!addMuon) { //No additional muons found
				h_e_tau_cutFlow->Fill("1 e & 0 #mu", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax
							&& tmpJet->Charge != tmpElectron->Charge) { //Quality  OS tau
						taus.push_back(i);
					}
				}
				if (taus.size() >= 1) {//Quality tau
					h_e_tau_cutFlow->Fill("Quality #tau", 1);
					v_tau_1 = tmpElectron->P4()
					tmpJet = (Jet*)branchJet->At(taus[0]);
					v_tau_0 = tmpJet->P4();
					tmpMPT = (MissingET*)branchMissingET->At(0);
					v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(branchParticle, branchElectron,
								branchMuon, branchJet,
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
					mPT_pT = tmpMPT->MET;
					mPT_phi = tmpMPT->Phi;
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
					tmpWeight = (Weight*)branchWeights->At(0);
					weight = tmpWeight->Weight;
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
		for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through Electons
			tmpElectron = (Electron*) branchElectron->At(i);
			if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
					&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electon
				addElectron = true;
				break;
			}
		}
		if (!addElectron) { //No additional electons
			for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
				tmpMuon = (Muon*) branchMuon->At(i);
				if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
						&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
					addMuon = true;
					break;
				}
			}
			if (!addMuon) { //No additional muons found
				h_tau_tau_b_b_cutFlow->Fill("0 e & 0 #mu", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax) { //Quality tau
						taus.push_back(i);
					}
				}
				if (taus.size() >= 2) {//2 quality taus
					h_tau_tau_b_b_cutFlow->Fill("Quality #tau#tau", 1);
					if (getOSTauTauPair(branchJet, &taus, &tau_0, &tau_1)) { //OS Tau pair
						tmpJet = (Jet*)branchJet->At(tau_0);
						v_tau_0 = tmpJet->P4();
						tmpJet = (Jet*)branchJet->At(tau_1);
						v_tau_1 = tmpElectron->P4();
						tmpMPT = (MissingET*)branchMissingET->At(0);
						v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
						gen_mctMatch = false;
						if (options["-t"] == "1") {
							gen_mctMatch = truthCut(branchParticle, branchElectron,
								branchMuon, branchJet,
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
						mPT_pT = tmpMPT->MET;
						mPT_phi = tmpMPT->Phi;
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
						tmpWeight = (Weight*)branchWeights->At(0);
						weight = tmpWeight->Weight;
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
		for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
			tmpMuon = (Muon*) branchMuon->At(i);
			if (tmpMuon->PT > muPTMin && std::abs(tmpMuon->Eta) < muEtaMax
					&& tmpMuon->IsolationVar < muIsoMax) { //Quality muon
				muons.push_back(i);
			} else if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
					&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 2 && !addMuon && 
			((Muon*)branchMuon->At(muons[0]))->Charge != ((Muon*)branchMuon->At(muons[1]))->Charge) { //One quality OS muon pair found and no additional muons
			h_mu_mu_cutFlow->Fill("Quality #mu#mu", 1);
			for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through electrons
				tmpElectron = (Electron*) branchElectron->At(i);
				if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
						&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electron
					addElectron = true;
					break;
				}
			}
			if (!addElectron) { //No additional electrons found
				h_mu_mu_cutFlow->Fill("2 #mu & 0 e", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax) { //Quality tau
						taus.push_back(i);
						break;
					}
				}
				if (taus.size() == 0) {//No taus
					h_mu_mu_cutFlow->Fill("0 #tau", 1);
					tmpMuon = (Muon*) branchMuon->At(muons[0]);
					v_tau_0 = tmpMuon->P4();
					tmpMuon = (Muon*) branchMuon->At(muons[1]);
					v_tau_1 = tmpMuon->P4();
					tmpMPT = (MissingET*)branchMissingET->At(0);
					v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(branchParticle, branchElectron,
								branchMuon, branchJet,
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
					mPT_pT = tmpMPT->MET;
					mPT_phi = tmpMPT->Phi;
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
					tmpWeight = (Weight*)branchWeights->At(0);
					weight = tmpWeight->Weight;
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
		for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
			tmpMuon = (Muon*) branchMuon->At(i);
			if (tmpMuon->PT > muPTMin && std::abs(tmpMuon->Eta) < muEtaMax
					&& tmpMuon->IsolationVar < muIsoMax) { //Quality muon
				muons.push_back(i);
			} else if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
					&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
				addMuon = true;
				break;
			}
		}
		if (muons.size() == 1 && !addMuon) { //One quality muon found and no additional muons
			h_e_mu_cutFlow->Fill("Quality #mu", 1);
			for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through Electons
				tmpElectron = (Electron*) branchElectron->At(i);
				if (tmpElectron->PT > ePTMin && std::abs(tmpElectron->Eta) < eEtaMax
						&& tmpElectron->IsolationVar < eIsoMax) { //Quality Electron
					electrons.push_back(i);
				} else if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
						&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electon
					addElectron = true;
					break;
				}
			}
			if (electrons.size() == 1 && !addElectron) { //No additional electrons found
				h_e_mu_cutFlow->Fill("1 #mu & 1 e", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax) { //Quality tau
						taus.push_back(i);
						break;
					}
				}
				if (taus.size() == 0) {//No taus
					h_e_mu_cutFlow->Fill("0 #tau", 1);
					tmpMuon = (Muon*) branchMuon->At(muons[0]);
					v_tau_0 = tmpMuon->P4();
					tmpElectron = (Electon*) branchElectron->At(electrons[1]);
					v_tau_1 = tmpElectron->P4();
					tmpMPT = (MissingET*)branchMissingET->At(0);
					v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(branchParticle, branchElectron,
								branchMuon, branchJet,
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
					mPT_pT = tmpMPT->MET;
					mPT_phi = tmpMPT->Phi;
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
					tmpWeight = (Weight*)branchWeights->At(0);
					weight = tmpWeight->Weight;
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
		for (int i = 0; i < branchElectron->GetEntries(); i++) { //Loop through Electons
			tmpElectron = (Electron*) branchElectron->At(i);
			if (tmpElectron->PT > ePTMin && std::abs(tmpElectron->Eta) < eEtaMax
					&& tmpElectron->IsolationVar < eIsoMax) { //Quality Electron
				electrons.push_back(i);
			} else if (tmpElectron->PT > ePTMinAdd && std::abs(tmpElectron->Eta) < eEtaMaxAdd
					&& tmpElectron->IsolationVar < eIsoMaxAdd) { //Additional electon
				addElectron = true;
				break;
			}
		}
		if (electrons.size() == 2 && !addElectron && 
			((Electron*)branchElectron->At(electrons[0]))->Charge != ((Electron*)branchElectron->At(electrons[1]))->Charge) { //One quality OS electron pair found and no additional electrons
			h_e_e_cutFlow->Fill("Quality ee", 1);
			for (int i = 0; i < branchMuon->GetEntries(); i++) { //Loop through muons
				tmpMuon = (Muon*) branchMuon->At(i);
				if (tmpMuon->PT > muPTMinAdd && std::abs(tmpMuon->Eta) < muEtaMaxAdd
						&& tmpMuon->IsolationVar < muIsoMaxAdd) { //Additional muon
					addMuon = true;
					break;
				}
			}
			if (!addMuon) { //No additional muons found
				h_e_e_cutFlow->Fill("2 e & 0 #mu", 1);
				for (int i = 0; i < branchJet->GetEntries(); i++) { //Loop through jets
					tmpJet = (Jet*) branchJet->At(i);
					if (tmpJet->TauTag == 1 && tmpJet->BTag == 0 && tmpJet->PT > tauPTMin
							&& std::abs(tmpJet->Eta) < tauEtaMax) { //Quality tau
						taus.push_back(i);
						break;
					}
				}
				if (taus.size() == 0) {//No taus
					h_e_e_cutFlow->Fill("0 #tau", 1);
					tmpElectron = (Electron*) branchElectron->At(electrons[0]);
					v_tau_0 = tmpElectron->P4();
					tmpElectron = (Electron*) branchElectron->At(electrons[1]);
					v_tau_1 = tmpElectron->P4();
					tmpMPT = (MissingET*)branchMissingET->At(0);
					v_higgs_tt = getHiggs2Taus(tmpMPT, v_tau_0, v_tau_1);
					gen_mctMatch = false;
					if (options["-t"] == "1") {
						gen_mctMatch = truthCut(branchParticle, branchElectron,
								branchMuon, branchJet,
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
					mPT_pT = tmpMPT->MET;
					mPT_phi = tmpMPT->Phi;
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
					tmpWeight = (Weight*)branchWeights->At(0);
					weight = tmpWeight->Weight;
					e_e->Fill();
					h_datasetSizes->Fill("e e", 1);
					eventAccepted = true;
				}
			}
		}
	}
	//___________________________________
	std::cout << "Event loop complete\n";
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