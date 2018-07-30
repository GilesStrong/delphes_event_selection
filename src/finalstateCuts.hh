/*
 * finalstateCuts.hh
 *
 *  Created on: 8 Apr 2016
 *      Author: giles
 */

#ifndef FINALSTATECUTS_HH_
#define FINALSTATECUTS_HH_

double ePTMin = 20.0; //Gev
double eEtaMax = 2.4;
double eIsoMax = 0.2;
double muPTMin = 20.0; //Gev
double muEtaMax = 2.4;
double muIsoMax = 0.2;
double ePTMinAdd = 20.0; //Gev
double eEtaMaxAdd = 2.4;
double eIsoMaxAdd = 0.2;
double muPTMinAdd = 20.0; //Gev
double muEtaMaxAdd = 2.4;
double muIsoMaxAdd = 0.2;
double tauPTMin; //Gev
double tauEtaMax;

bool massCut = false;

void finalstateSet(std::string set) {
	if (set == "e_tau_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.4;
	} else if (set == "mu_tau_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.4;
	} else if (set == "tau_tau_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.4;
	} else if (set == "e_e_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.4;
	} else if (set == "e_mu_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.4;
	} else if (set == "mu_mu_b_b") {
		tauPTMin = 20.0;
		tauEtaMax = 2.4;
	}
}
#endif /* FINALSTATECUTS_HH_ */
