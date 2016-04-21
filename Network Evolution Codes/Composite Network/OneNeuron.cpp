/*
 * OneNeuron.cpp
 *
 *  Created on: 2010/05/27
 *      Author: Amit Lal & Yoshitaka Oku
 */

#include "OneNeuron.h"
#include "math.h"
#include <ctime>
#include <cstdlib>
//#include "matrix.h"

namespace pBC {

OneNeuron::OneNeuron() {
	// TODO Auto-generated constructor stub
	Initialize();
	//cout << "initialization was done!\n";
}

OneNeuron::~OneNeuron() {
	// TODO Auto-generated destructor stub
}

/*
array &operator +(array x, array y) {
   if(y.size() != x.size()) {
      cerr << "array operator *(array, array):Invalid array size.\n";
      exit(EXIT_FAILURE);
   }
   static array z(x.size());
   for(int i = 0;i < x.size();i++) {
	   z[i] = x[i]+y[i];
   }
   return z;
}
*/

void OneNeuron::Initialize()
{
	// parameter initialization
	C = 21; //pF
	gNa = 28; //ns
	gK = 11.2; //ns

	ENa = 50; //mV
	EK = -85; //mV
	EL = -70; //mV
	Esyn_e = 0; //mV
	ECL = -90; //mV

	theta_mNa = -34; //mV
	sigma_mNa = -5; //mV
	theta_n = -29; //mV
	sigma_n = -4; //mV
	Tau_n = 10; //ms
	theta_h = -53; //mV
	sigma_h = 6; //mV
	Tau_h = 10000; //ms
	theta_s = -10; //mV
	sigma_s = -5; //mV
	Tau_s = 5; //ms
	theta_mNaP = -45.1; //mV
	sigma_mNaP = -5; //mV
}

void OneNeuron::GNaPo_set( double value ) {
	gNaP = value;
}

void OneNeuron::GLeako_set( double value ) {
	gL = value;
}

double OneNeuron::GNaP_return()
{
	return gNaP;
}

double OneNeuron::EL_return()
{
	return EL;
}

int OneNeuron::RN1_return()
{
	return RN1;
}

int OneNeuron::RN2_return()
{
	return RN2;
}

void OneNeuron::Gtonic_e(double ggtonic_e)
{
	gtonic_e = ggtonic_e;
}

double OneNeuron::Vrate(double Vv, double nv, double hv,
											double sum_gsyn_eXsv, double sum_gsyn_iXsv)
{
	double I_Na, I_K, I_NaP, I_L, I_tonic_e, I_syn, vDot;
	m_inf_Na = 1/(1+exp((Vv-theta_mNa)/sigma_mNa));
	m_inf_NaP = 1/(1+exp((Vv-theta_mNaP)/sigma_mNaP));

	I_Na = gNa*pow(m_inf_Na,3)*(1-nv)*(Vv-ENa);
	I_K = gK*pow(nv,4)*(Vv-EK);
	I_NaP = gNaP*m_inf_NaP*hv*(Vv-ENa);
	I_L = gL*(Vv-EL);
	I_tonic_e = gtonic_e*(Vv-Esyn_e);
	I_syn = sum_gsyn_eXsv*(Vv-Esyn_e) +  sum_gsyn_iXsv*(Vv-ECL);

	vDot = 1/C*( -I_NaP - I_Na - I_K - I_L - I_tonic_e - I_syn);
	return vDot;
}

double OneNeuron::nrate(double Vn, double nn)
{
	double n_inf, tau_n, nDot;
	n_inf = 1/(1+exp((Vn-theta_n)/sigma_n));
	tau_n = Tau_n/cosh((Vn-theta_n)/(2*sigma_n));
	nDot = (n_inf - nn)/tau_n;
	return nDot;
}

double OneNeuron::hrate(double Vh, double hh)
{
	double h_inf, tau_h, hDot;
	h_inf = 1/(1+exp((Vh-theta_h)/sigma_h));
	tau_h = Tau_h/cosh((Vh-theta_h)/(2*sigma_h));
	hDot = (h_inf - hh)/tau_h;
	return hDot;
}

double OneNeuron::srate(double Vs, double ss)
{
	double s_inf, sDot;
	s_inf = 1/(1+exp((Vs-theta_s)/sigma_s));
	sDot = ( (1-ss)*s_inf - ss )/Tau_s;
	return sDot;
}

double OneNeuron::SF(double V, double n, double h,
		double sum_gsyn_eXsv, double ddv)
{
	double m, mt, mdv, mm, mmt, mmdv;
	double J111, J11;

	m = 1/(1+exp((V-theta_mNa)/sigma_mNa));
	mt = 1/(1+exp((V+ddv-theta_mNa)/sigma_mNa));
	mdv = (mt - m) / ddv;

	mm = 1/(1+exp((V-theta_mNaP)/sigma_mNaP));
	mmt = 1/(1+exp((V+ddv-theta_mNaP)/sigma_mNaP));
	mmdv = (mmt - mm) / ddv;

    J111 = gNa*pow(m,3)*(1-n)+gNaP*mm*h
    		+(gNa*3*pow(m,2)*mdv*(1-n)+gNaP*mmdv*h)*(V-ENa);
    J11 = -(J111 + gK*pow(n,4) + gL + sum_gsyn_eXsv)/C;
    return J11;
}

/*
void OneNeuron::Euler_LL
	(double V_now, double n_now, double h_now,
			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step)
{
    double ddv = 1; // mesh of v at numerical differentiation
    double dtmin = 0.1;  // time mesh at unstable state (msec)
    double thresh = 0.0001; // Threshold for judging stable/unstable state
     double thresh_ext=0.07;
	double V, n, h, s, sum_gsyn_eXs, sum_gsyn_iXs, dt;
	double m, mm, taun, tauh;
	double n8, n80, taun00, n8dv, taundv;
	double h8, h80, tauh00, h8dv, tauhdv;
	//double s8;
	double J11n;
	int ndivide = round(t_step / dtmin);
	const int size = 3;
	MAT::matrix FXDX(size,size), ATF1(size,size);
	array x(size), xt(size), xdot(size);

	V = V_now;
	n = n_now;
	h = h_now;
	s = s_now;
	sum_gsyn_eXs = ext_gsyn_eXs;
	sum_gsyn_iXs = ext_gsyn_iXs;
	dt = t_step;

	double V_dot = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
	double n_dot = nrate(V, n);
	double h_dot = hrate(V, h);
	double s_dot = srate(V, s);

	// Euler Method
	double V1 = V + V_dot*dt;
	double n1 = n + n_dot*dt;
	double h1 = h + h_dot*dt;
	double s1 = s + s_dot*dt;

	// Judging the stability
	J11n = fabs(SF(V1,n1,h1,sum_gsyn_eXs,ddv)-SF(V,n,h,sum_gsyn_eXs,ddv));
	if ((J11n > thresh)|| (sum_gsyn_eXs > thresh_ext)){
		// recomputing
	    dt = dtmin;
		x[0]=V;
		x[1]=n;
		x[2]=h;
		for (int jjj=0; jjj<ndivide; jjj++){
			xt = x;
			m = 1/(1+exp((V-theta_mNa)/sigma_mNa));
			mm = 1/(1+exp((V-theta_mNaP)/sigma_mNaP));
			n8 = 1/(1+exp((V-theta_n)/sigma_n));
			taun = Tau_n/cosh((V-theta_n)/(2*sigma_n));
			h8 = 1/(1+exp((V-theta_h)/sigma_h));
			tauh = Tau_h/cosh((V-theta_h)/(2*sigma_h));
			//s8 = 1/(1+exp((V-theta_s)/sigma_s));

			xdot[0] = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
			xdot[1] = nrate(V, n);
			xdot[2] = hrate(V, h);

			// Local linearization
			n80 = 1/(1+exp((V+ddv-theta_n)/sigma_n));
		    taun00 = Tau_n / cosh((V + ddv -theta_n)/(2*sigma_n));
		    n8dv = (n80/taun00-n8/taun)/ddv;
		    taundv = (1/taun00-1/taun)/ddv;

			h80 = 1/(1+exp((V+ddv-theta_h)/sigma_h));
			tauh00 = Tau_h/cosh((V+ddv-theta_h)/(2*sigma_h));
		    h8dv = (h80/tauh00-h8/tauh)/ddv;
		    tauhdv = (1/tauh00-1/tauh)/ddv;

			// start Jacobian
		    // F/v
		    FXDX(0,0) = SF(V+ddv, n, h, sum_gsyn_eXs, ddv);
		    FXDX(0,1) = (gNa*pow(m,3)*(V-ENa)-4*gK*pow(n,3)*(V-EK))/C;
		    FXDX(0,2) = -gNaP*mm*(V-ENa)/C;
			//F/n
		    FXDX(1,0) = n8dv -n*taundv;
		    FXDX(1,1) = -1/taun;
		    FXDX(1,2) = 0;
			//F/h
		    FXDX(2,0) = h8dv -h*tauhdv;
		    FXDX(2,1) = 0;
		    FXDX(2,2) = -1/tauh;

		    // inv(FXDX)*(expm(FXDX*dt)-I)=Idt+FXDX*dt^2/2+  +FXDX^N-1/N!
		    MAT::matrix AA(size,size), AAA(size,size), IU3(size,size);
		    AA=FXDX*dt;
		    IU3.identity();
		    AAA=IU3*dt;
		    ATF1=AAA;
		    int LL=1;
		    double maxAAA;
		    for (int n1=1; n1<=5; n1++) {
		    	maxAAA=max(AAA);
				if (maxAAA>0.0000001) {
					for (int n2=1; n2<=8; n2++) {
						LL=LL+1;
						AAA=AAA*AA/(double)LL;
						ATF1=ATF1+AAA;
					}
				}
		    }
		    x=(ATF1*xdot)+xt;
			V = x[0];
			n = x[1];
			h = x[2];
		}
		V_next = V;
		n_next = n;
		h_next = h;
		s_next = s1;
	}
	else {
		V_next = V1;
		n_next = n1;
		h_next = h1;
		s_next = s1;
	}
}
*/

void OneNeuron::Euler_RK4
	(double V_now, double n_now, double h_now,
			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step)
{
    double ddv = 1; // mesh of v at numerical differentiation
    double dtmin = 0.05;  // time mesh at unstable state (msec)
    double thresh = 0.0001; // Threshold for judging stable/unstable state
    double thresh_ext=0.07;
	double V, n, h, s, sum_gsyn_eXs, sum_gsyn_iXs, dt;
	double J11n;
	int ndivide = round(t_step / dtmin);

	V = V_now;
	n = n_now;
	h = h_now;
	s = s_now;
	sum_gsyn_eXs = ext_gsyn_eXs;
	sum_gsyn_iXs = ext_gsyn_iXs;
	dt = t_step;

	double V_dot = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
	double n_dot = nrate(V, n);
	double h_dot = hrate(V, h);
	double s_dot = srate(V, s);

	// Euler Method
	double V1 = V + V_dot*dt;
	double n1 = n + n_dot*dt;
	double h1 = h + h_dot*dt;
	double s1 = s + s_dot*dt;

	// Judging the stability
	J11n = fabs(SF(V1,n1,h1,sum_gsyn_eXs,ddv)-SF(V,n,h,sum_gsyn_eXs,ddv));
	if ((J11n > thresh)|| (sum_gsyn_eXs > thresh_ext)){
		// recomputing
		dt = dtmin;
		for (int jjj=0; jjj<ndivide; jjj++){
			// Euler Method
			//V_dot = Vrate(V, n, h, sum_gsyn_eXs);
			//n_dot = nrate(V, n);
			//h_dot = hrate(V, h);
			//h_dot = srate(V, s);

			//V += V_dot*dt;
			//n += n_dot*dt;
			//h += h_dot*dt;
			//s += s_dot*dt;

			// Runge-Kutta Method
			double V0_dot = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
			double n0_dot = nrate(V, n);
			double h0_dot = hrate(V, h);
			double s0_dot = srate(V, s);

			V1 = V + V0_dot*dt/2.0;
			n1 = n + n0_dot*dt/2.0;
			h1 = h + h0_dot*dt/2.0;
			s1 = s + s0_dot*dt/2.0;

			double V1_dot = Vrate(V1, n1, h1, sum_gsyn_eXs, sum_gsyn_iXs);
			double n1_dot = nrate(V1, n1);
			double h1_dot = hrate(V1, h1);
			double s1_dot = srate(V1, s1);

			double V12 = V + V1_dot*dt/2.0;
			double n12 = n + n1_dot*dt/2.0;
			double h12 = h + h1_dot*dt/2.0;
			double s12 = s + s1_dot*dt/2.0;

			double V12_dot = Vrate(V12, n12, h12, sum_gsyn_eXs, sum_gsyn_iXs);
			double n12_dot = nrate(V12, n12);
			double h12_dot = hrate(V12, h12);
			double s12_dot = srate(V12, s12);

			double V2 = V + V12_dot*dt;
			double n2 = n + n12_dot*dt;
			double h2 = h + h12_dot*dt;
			double s2 = s + s12_dot*dt;

			double V2_dot = Vrate(V2, n2, h2, sum_gsyn_eXs, sum_gsyn_iXs);
			double n2_dot = nrate(V2, n2);
			double h2_dot = hrate(V2, h2);
			double s2_dot = srate(V2, s2);

			V += dt/6*( V0_dot + 2*V1_dot + 2*V12_dot + V2_dot);
			n += dt/6*( n0_dot + 2*n1_dot + 2*n12_dot + n2_dot);
			h += dt/6*( h0_dot + 2*h1_dot + 2*h12_dot + h2_dot);
			s += dt/6*( s0_dot + 2*s1_dot + 2*s12_dot + s2_dot);
		}
		V_next = V;
		n_next = n;
		h_next = h;
		s_next = s;
	}
	else {
		V_next = V1;
		n_next = n1;
		h_next = h1;
		s_next = s1;
	}
}

void OneNeuron::RK4_RK4
	(double V_now, double n_now, double h_now,
			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step)
{
    double ddv = 1; // mesh of v at numerical differentiation
    double dtmin = 0.05;  // time mesh at unstable state (msec)
    double thresh = 0.0001; // Threshold for judging stable/unstable state
    double thresh_ext=0.07;
	double V, n, h, s, sum_gsyn_eXs, sum_gsyn_iXs, dt;
	double J11n;
	int ndivide = round(t_step / dtmin);

	V = V_now;
	n = n_now;
	h = h_now;
	s = s_now;
	sum_gsyn_eXs = ext_gsyn_eXs;
	sum_gsyn_iXs = ext_gsyn_iXs;
	dt = t_step;

	// Runge-Kutta Method
	double V0_dot = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
	double n0_dot = nrate(V, n);
	double h0_dot = hrate(V, h);
	double s0_dot = srate(V, s);

	double V1 = V + V0_dot*dt/2.0;
	double n1 = n + n0_dot*dt/2.0;
	double h1 = h + h0_dot*dt/2.0;
	double s1 = s + s0_dot*dt/2.0;

	double V1_dot = Vrate(V1, n1, h1, sum_gsyn_eXs, sum_gsyn_iXs);
	double n1_dot = nrate(V1, n1);
	double h1_dot = hrate(V1, h1);
	double s1_dot = srate(V1, s1);

	double V12 = V + V1_dot*dt/2.0;
	double n12 = n + n1_dot*dt/2.0;
	double h12 = h + h1_dot*dt/2.0;
	double s12 = s + s1_dot*dt/2.0;

	double V12_dot = Vrate(V12, n12, h12, sum_gsyn_eXs, sum_gsyn_iXs);
	double n12_dot = nrate(V12, n12);
	double h12_dot = hrate(V12, h12);
	double s12_dot = srate(V12, s12);

	double V2 = V + V12_dot*dt;
	double n2 = n + n12_dot*dt;
	double h2 = h + h12_dot*dt;
	double s2 = s + s12_dot*dt;

	double V2_dot = Vrate(V2, n2, h2, sum_gsyn_eXs, sum_gsyn_iXs);
	double n2_dot = nrate(V2, n2);
	double h2_dot = hrate(V2, h2);
	double s2_dot = srate(V2, s2);

	V1 = V + dt/6*( V0_dot + 2*V1_dot + 2*V12_dot + V2_dot);
	n1 = n + dt/6*( n0_dot + 2*n1_dot + 2*n12_dot + n2_dot);
	h1 = h + dt/6*( h0_dot + 2*h1_dot + 2*h12_dot + h2_dot);
	s1 = s + dt/6*( s0_dot + 2*s1_dot + 2*s12_dot + s2_dot);

	// Judging the stability
	J11n = fabs(SF(V1,n1,h1,sum_gsyn_eXs,ddv)-SF(V,n,h,sum_gsyn_eXs,ddv));
	if ((J11n > thresh)|| (sum_gsyn_eXs > thresh_ext)){
		// recomputing
		dt = dtmin;
		for (int jjj=0; jjj<ndivide; jjj++){
			// Runge-Kutta Method
			V0_dot = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
			n0_dot = nrate(V, n);
			h0_dot = hrate(V, h);
			s0_dot = srate(V, s);

			V1 = V + V0_dot*dt/2.0;
			n1 = n + n0_dot*dt/2.0;
			h1 = h + h0_dot*dt/2.0;
			s1 = s + s0_dot*dt/2.0;

			V1_dot = Vrate(V1, n1, h1, sum_gsyn_eXs, sum_gsyn_iXs);
			n1_dot = nrate(V1, n1);
			h1_dot = hrate(V1, h1);
			s1_dot = srate(V1, s1);

			V12 = V + V1_dot*dt/2.0;
			n12 = n + n1_dot*dt/2.0;
			h12 = h + h1_dot*dt/2.0;
			s12 = s + s1_dot*dt/2.0;

			V12_dot = Vrate(V12, n12, h12, sum_gsyn_eXs, sum_gsyn_iXs);
			n12_dot = nrate(V12, n12);
			h12_dot = hrate(V12, h12);
			s12_dot = srate(V12, s12);

			V2 = V + V12_dot*dt;
			n2 = n + n12_dot*dt;
			h2 = h + h12_dot*dt;
			s2 = s + s12_dot*dt;

			V2_dot = Vrate(V2, n2, h2, sum_gsyn_eXs, sum_gsyn_iXs);
			n2_dot = nrate(V2, n2);
			h2_dot = hrate(V2, h2);
			s2_dot = srate(V2, s2);

			V += dt/6*( V0_dot + 2*V1_dot + 2*V12_dot + V2_dot);
			n += dt/6*( n0_dot + 2*n1_dot + 2*n12_dot + n2_dot);
			h += dt/6*( h0_dot + 2*h1_dot + 2*h12_dot + h2_dot);
			s += dt/6*( s0_dot + 2*s1_dot + 2*s12_dot + s2_dot);
		}
		V_next = V;
		n_next = n;
		h_next = h;
		s_next = s;
	}
	else {
		V_next = V1;
		n_next = n1;
		h_next = h1;
		s_next = s1;
	}
}

void OneNeuron::RK4
	(double V_now, double n_now, double h_now,
			double s_now, double ext_gsyn_eXs, double ext_gsyn_iXs, double t_step)
{
	double V, n, h, s, sum_gsyn_eXs, sum_gsyn_iXs, dt;
	V = V_now;
	n = n_now;
	h = h_now;
	s = s_now;
	sum_gsyn_eXs = ext_gsyn_eXs;
	sum_gsyn_iXs = ext_gsyn_iXs;
	dt = t_step;

	double V0_dot = Vrate(V, n, h, sum_gsyn_eXs, sum_gsyn_iXs);
	double n0_dot = nrate(V, n);
	double h0_dot = hrate(V, h);
	double s0_dot = srate(V, s);

	double V1 = V + V0_dot*dt/2.0;
	double n1 = n + n0_dot*dt/2.0;
	double h1 = h + h0_dot*dt/2.0;
	double s1 = s + s0_dot*dt/2.0;

	double V1_dot = Vrate(V1, n1, h1, sum_gsyn_eXs, sum_gsyn_iXs);
	double n1_dot = nrate(V1, n1);
	double h1_dot = hrate(V1, h1);
	double s1_dot = srate(V1, s1);

	double V12 = V + V1_dot*dt/2.0;
	double n12 = n + n1_dot*dt/2.0;
	double h12 = h + h1_dot*dt/2.0;
	double s12 = s + s1_dot*dt/2.0;

	double V12_dot = Vrate(V12, n12, h12, sum_gsyn_eXs, sum_gsyn_iXs);
	double n12_dot = nrate(V12, n12);
	double h12_dot = hrate(V12, h12);
	double s12_dot = srate(V12, s12);

	double V2 = V + V12_dot*dt;
	double n2 = n + n12_dot*dt;
	double h2 = h + h12_dot*dt;
	double s2 = s + s12_dot*dt;

	double V2_dot = Vrate(V2, n2, h2, sum_gsyn_eXs, sum_gsyn_iXs);
	double n2_dot = nrate(V2, n2);
	double h2_dot = hrate(V2, h2);
	double s2_dot = srate(V2, s2);

	V_next = V + dt/6*( V0_dot + 2*V1_dot + 2*V12_dot + V2_dot);
	n_next = n + dt/6*( n0_dot + 2*n1_dot + 2*n12_dot + n2_dot);
	h_next = h + dt/6*( h0_dot + 2*h1_dot + 2*h12_dot + h2_dot);
	s_next = s + dt/6*( s0_dot + 2*s1_dot + 2*s12_dot + s2_dot);
}

double OneNeuron::V_Next()
{
	return V_next;
}

double OneNeuron::n_Next()
{
	return n_next;
}

double OneNeuron::h_Next()
{
	return h_next;
}

double OneNeuron::s_Next()
{
	return s_next;
}

}
