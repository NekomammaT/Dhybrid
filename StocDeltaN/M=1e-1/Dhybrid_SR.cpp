#include "source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "Dhybrid_SR" // model name

// ---------- for PDE ----------------------
#define MAXSTEP 100000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define AS 2.1e-9
#define MM 0.1
#define PHIC (0.1*sqrt(2))
#define MU2 10 //10.7
// -----------------------------------------

// ---------- box size & step h ------------
#define PSIMIN -(1e-3)
#define PSIMAX (1e-3)
#define HPSIOPSI (1e-2) // hpsi/|psi|
// -----------------------------------------

// ---------- for SDE ----------------------
#define RECURSION 100000 // recursion for power spectrum
#define TIMESTEP (1e-3) // time step : delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.5 // 0.1 // calc. PS every DELTAN e-folds
// -----------------------------------------


int main(int argc, char** argv)
{
  // -------- for Dhybrid ------------
  if (argc != 4) {
    cout << "Pi2, Dwater, DeltaVend/V0を正しく指定してください" << endl;
    return 1;
  }
  
  double Pi2 = atof(argv[1]);
  double Dwater = atof(argv[2]);
  double DeltaVoV0 = atof(argv[3]);
  
  double mu1 = Pi2/MM/MM/PHIC;
  double Lambda4 = AS*12*M_PI*M_PI/mu1/mu1;
  double sigmapsi = sqrt(Dwater*Lambda4*sqrt(Pi2)/48./sqrt(2*M_PI)/M_PI);
  
  double phimin = 0.1409; //0.141; // - 1./mu1/10;

  double hphi = 1./mu1/10;
  double hpsimin = sigmapsi/10;
  double psiin = sigmapsi;
  double phimax, phiin;

  phimax = PHIC + 20./mu1;
  phiin = PHIC + 15./mu1;

  double Nmax = 10;
  // ---------------------------------
  
  // ---------- start stop watch ----------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // --------------------------------------
  

  // ---------- set box ---------------
  double h = hphi, //HPHI,
    sitev = phimin; //PHIMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= //PHIMAX
	 phimax) {
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();
  
  sitev = PSIMIN;
  while (sitev <= PSIMAX) {
    h = max(fabs(sitev)*HPSIOPSI, hpsimin //HPSIMIN
	    );

    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  site.clear();

  sitepack.push_back(xsite);
  xsite.clear();
  // ----------------------------------

  vector<double> params;
  double rhoc;
  vector< vector<double> > xpi = {{phiin, //PHIIN,
      psiin //PSIIN
    }}; // set i.c. for inflationary trajectories

  string model = string(MODEL) + string("_Pi2=") + to_string((int)Pi2)
    + string("_D=") + to_string((int)Dwater);

  rhoc = (1-DeltaVoV0) * Lambda4;

  params = {MAXSTEP,TOL,2,rhoc
	    ,(double)sitepack[0].size(),TIMESTEP,Nmax //NMAX
	    ,DELTAN,RECURSION //,NCUT
	    ,Pi2,Dwater,Lambda4,mu1,0 // for Dhybrid
  }; // set parameters
  
  StocDeltaN sdn(model,sitepack,xpi,0,params); // declare the system
  //sdn.solve_fg();

  sdn.sample();
  cout << endl;

  Nmax = sdn.return_t() - 5;
  params[6] = Nmax;

  StocDeltaN sdn2(model,sitepack,xpi,0,params);
  string Mnfile = string("./data/Mn_") + string(MODEL) + string ("_Pi2=") + to_string((int)Pi2)
  + string("_D=") + to_string((int)Dwater) + string(".dat");
  sdn2.import_fg(Mnfile);
  sdn2.solve();

  
  // ---------- stop stop watch ----------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------

  cout << endl;
}


// ---------- Lagrangian params. X[0]=phi, X[1]=psi ----------

double StocDeltaN::V(vector<double> &X)
{
  return //LAMBDA4
    Lambda4*((1-X[1]*X[1]/MM/MM)*(1-X[1]*X[1]/MM/MM)
	     + 2*X[0]*X[0]*X[1]*X[1]/PHIC/PHIC/MM/MM + (X[0]-PHIC)/mu1 //MU1
	     - (X[0]-PHIC)*(X[0]-PHIC)/MU2/MU2);
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  if (I == 0) {
    return //LAMBDA4
      Lambda4*(1./mu1 //MU1
	       - 2*(X[0]-PHIC)/MU2/MU2 + 4*X[0]*X[1]*X[1]/PHIC/PHIC/MM/MM);
  } else {
    return //LAMBDA4
      Lambda4*(4*X[0]*X[0]*X[1]/PHIC/PHIC/MM/MM
	       - 4*X[1]*(1-X[1]*X[1]/MM/MM)/MM/MM);
  }
}

double StocDeltaN::VIJ(vector<double> &X, int I, int J)
{
  if (I == 0 && J == 0) {
    return //LAMBDA4
      Lambda4*(4*X[1]*X[1]/MM/MM/PHIC/PHIC - 2./MU2/MU2);
  } else if (I == 1 && J == 1) {
    return 4* //LAMBDA4
      Lambda4*(MM*MM*(X[0]*X[0]-PHIC*PHIC) + 3*PHIC*PHIC*X[1]*X[1])
      /MM/MM/MM/MM/PHIC/PHIC;
  } else {
    return 8*X[0]*X[1]* //LAMBDA4
      Lambda4/MM/MM/PHIC/PHIC;
  }
}

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

/*
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
      + V(psv[0])*(Dwater-1)/abs(psv[0][I])/24./M_PI/M_PI;
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


bool StocDeltaN::EndSurface(vector< vector<double> > &psv)
{
  if (SRend) {
    return abs(etaV(psv[0])) < 2 || psv[0][0] > PHIC;
    //return eV(psv[0]) < 1e-5*(50./PI2)*(50./PI2);
  } else {
    return V(psv[0]) >= rhoc;
  }
}
