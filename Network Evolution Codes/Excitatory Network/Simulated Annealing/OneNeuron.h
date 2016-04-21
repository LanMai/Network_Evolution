/*
 * OneNeuron.h
 *
 *  Created on: 2010/05/27
 *      Author: Amit Lal & Yoshitaka Oku
 */

#ifndef ONENEURON_H_
#define ONENEURON_H_

namespace pBC {

class OneNeuron {

	double C; //pF
	double gNa; //ns
	double gK; //ns
	double gL; //ns

	double ENa; //mV
	double EK; //mV
	double ELo; //mV
	double Esyn_e; //mV
	double ECL; //mV

	double theta_mNa; //mV
	double sigma_mNa; //mV
	double theta_n; //mV
	double sigma_n; //mV
	double Tau_n; //ms
	double theta_h; //mV
	double sigma_h; //mV
	double Tau_h; //ms
	double theta_s; //mV
	double sigma_s; //mV
	double Tau_s; //ms
	double theta_mNaP; //mV
	double sigma_mNaP; //mV

	double m_inf_Na;
	double m_inf_NaP;

	double gtonic_e;
	double gNaPo, gNaP, EL;

	double V_next, n_next, h_next, s_next;

	int RN1;
	int RN2;

public:
	OneNeuron();
	virtual ~OneNeuron();

	void Initialize();
	void GNaPo_set( double value );
	void GLeako_set( double value );
	double GNaP_return();
	double EL_return();
	int RN1_return();
	int RN2_return();
	void Gtonic_e(double ggtonic_e);
	double Vrate(double Vv, double nv, double hv,
			double sum_gsyn_eXsv, double sum_gsyn_iXsv);
	double nrate(double Vn, double nn);
	double hrate(double Vh, double hh);
	double srate(double Vs, double ss);
	double SF(double V, double n, double h, double sum_gsyn_eXsv, double ddv);
//	void Euler_LL(double V_now, double n_now, double h_now,
//			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step);
	void Euler_RK4(double V_now, double n_now, double h_now,
			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step);
	void RK4_RK4(double V_now, double n_now, double h_now,
				double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step);
	void RK4(double V_now, double n_now, double h_now,
			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step);
	double V_Next();
	double n_Next();
	double h_Next();
	double s_Next();
};

}

#endif /* ONENEURON_H_ */
