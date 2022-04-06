#include "source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "Dhybrid_SR" // model name

// ---------- box size & step h ------------
#define PHIMIN 0.1409
#define PHIMAX 0.142
#define PSIMIN -(1e-3)
#define PSIMAX (1e-3)
#define HPHI (1e-5)
#define HPSIOPSI (1e-2) // hpsi/|psi|
#define HPSIMIN (1e-10)
//#define NCUT 50

#define DWATER 1000 // # of waterfalls
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define AS (2.189e-9)
#define PI2 50
#define MM 0.1
#define PHIC (0.1*sqrt(2))
#define MU2 11
#define MU1 (PI2/MM/MM/PHIC)
#define LAMBDA4 (AS*12*M_PI*M_PI/MU1/MU1)
// -----------------------------------------

#define RHOC (2.074038e-16) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 10000 // recursion for power spectrum
#define PHIIN 0.1418 // i.c. for phi
#define PSIIN 1e-10 //0 // i.c. for psi
#define TIMESTEP (1e-3) // time step : delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.5 // 0.1 // calc. PS every DELTAN e-folds
#define NMAX 25 // 28 // calc. PS for 0--NMAX e-folds
// -----------------------------------------


int main(int argc, char** argv)
{
  if (argc != 2) {
    cout << "Dwater を正しく指定してください" << endl;
    return 1;
  }
  
  // ---------- start stop watch ----------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // --------------------------------------


  // ---------- set box ---------------
  double h = HPHI, sitev = PHIMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= PHIMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitev = PSIMIN;
  while (sitev <= PSIMAX) {
    h = max(fabs(sitev)*HPSIOPSI,HPSIMIN);

    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();
  // ----------------------------------

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,RECURSION //,NCUT
  }; // set parameters

  vector< vector<double> > xpi = {{PHIIN,PSIIN}}; // set i.c. for inflationary trajectories

  string model = string(MODEL) + string("_D=") + to_string(DWATER);
  StocDeltaN sdn(model,sitepack,xpi,0,params); // declare the system
  
  //sdn.sample(); // obtain 1 sample path
  //sdn.sample_logplot(); // plot obtained sample path
  
  sdn.solve(); // solve PDE & SDE to obtain power spectrum
  sdn.f_logplot(0); // show plot of <N>
  sdn.f_logplot(1); // show plot of <delta N^2>
  sdn.calP_plot(); // show plot of power spectrum of zeta


  // ---------- stop stop watch ----------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}


// ---------- Lagrangian params. X[0]=phi, X[1]=psi ----------

double StocDeltaN::V(vector<double> &X)
{
  return LAMBDA4*((1-X[1]*X[1]/MM/MM)*(1-X[1]*X[1]/MM/MM)
		  + 2*X[0]*X[0]*X[1]*X[1]/PHIC/PHIC/MM/MM + (X[0]-PHIC)/MU1
		  - (X[0]-PHIC)*(X[0]-PHIC)/MU2/MU2);
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  if (I == 0) {
    return LAMBDA4*(1./MU1 - 2*(X[0]-PHIC)/MU2/MU2 + 4*X[0]*X[1]*X[1]/PHIC/PHIC/MM/MM);
  } else {
    return LAMBDA4*(4*X[0]*X[0]*X[1]/PHIC/PHIC/MM/MM
		    - 4*X[1]*(1-X[1]*X[1]/MM/MM)/MM/MM);
  }
}

/*
double StocDeltaN::metric(vector<double> &X, int I, int J) // G_IJ
{
  if (I == J) {
    return 1;
  } else {
    return 0;
  }
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J) // G^IJ
{
  return metric(X,I,J);
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K) // \Gamma^I_JK
{
  return 0;
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L) // \partial_L \Gamma^I_JK
{
  return 0;
}
*/


double StocDeltaN::DI(int xp, int I, vector< vector<double> > &psv)
{
  if (I == 0) {
    return -VI(psv[0],I)/V(psv[0]);
  } else {
    return -VI(psv[0],I)/V(psv[0])
      + V(psv[0])*(DWATER-1)/abs(psv[0][I])/24./M_PI/M_PI;
  }
}

double StocDeltaN::DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv)
{
  if (I == J) {
    return V(psv[0])/12./M_PI/M_PI;
  } else {
    return 0;
  }
}

double StocDeltaN::gIa(int xp, int I, int alpha, vector< vector<double> > &psv)
{
  if (I == alpha) {
    return sqrt(V(psv[0])/3.)/2./M_PI;
  } else {
    return 0;
  }
}
