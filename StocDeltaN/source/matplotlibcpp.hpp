#ifndef INCLUDED_matplotlibcpp_hpp_
#define INCLUDED_matplotlibcpp_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

using namespace std;

class matplotlibcpp {
  FILE *p;

public:
  matplotlibcpp(){}
  void open();
  void xlabel(string label);
  void ylabel(string label);
  void yerrorbar(vector<double> X, vector<double> Y, vector<double> Yerr, double capsize, string fmt, double markersize, string ecolor, string markeredgecolor, string color);
  void plot(vector<double> X, vector<double> Y, double lw, string color);
  void xlog();
  void ylog();
  void contourf(vector<double> X, vector<double> Y, vector<double> Z, string label);
  void log_contourf(vector<double> X, vector<double> Y, vector<double> Z, string label);
  void show();
  void save(string filename);
  void close();
};

void matplotlibcpp::open()
{
  p = popen("python -c 'import code; import os; import sys; sys.stdout = sys.stderr = open(os.devnull, \"w\"); code.InteractiveConsole().interact()'", "w");
  fprintf(p, "import matplotlib.pyplot as plt\n");
}

void matplotlibcpp::xlabel(string label)
{
  fprintf(p, "plt.xlabel(\"%s\")\n", label.c_str());
}

void matplotlibcpp::ylabel(string label)
{
  fprintf(p, "plt.ylabel(\"%s\")\n", label.c_str());
}

void matplotlibcpp::yerrorbar(vector<double> X, vector<double> Y, vector<double> Yerr, double capsize, string fmt, double markersize, string ecolor, string markeredgecolor, string color)
{
  if (X.size() == Y.size() && Y.size() == Yerr.size()) {
    fprintf(p, "plt.errorbar([");

    bool first = true;
    for (int i=0; i<X.size(); i++) {
      if (isfinite(X[i]) && isfinite(Y[i]) && isfinite(Yerr[i])) {
	if (first) {
	  fprintf(p,"%e",X[i]);
	  first = false;
	} else {
	  fprintf(p, ",%e", X[i]);
	}
      }
    }

    fprintf(p, "], [");

    first = true;
    for (int i=0; i<Y.size(); i++) {
      if (isfinite(X[i]) && isfinite(Y[i]) && isfinite(Yerr[i])) {
	if (first) {
	  fprintf(p,"%e", Y[i]);
	  first = false;
	} else {
	  fprintf(p, ",%e", Y[i]);
	}
      }
    }

    fprintf(p, "], yerr = [");

    first = true;
    for (int i=0; i<Yerr.size(); i++) {
      if (isfinite(X[i]) && isfinite(Y[i]) && isfinite(Yerr[i])) {
	if (first) {
	  fprintf(p,"%e",Yerr[i]);
	  first = false;
	} else {
	  fprintf(p, ",%e", Yerr[i]);
	}
      }
    }
    
    fprintf(p, "], capsize=%f, fmt=\"%s\", markersize=%f, ecolor=\"%s\", markeredgecolor=\"%s\", color=\"%s\")\n", capsize, fmt.c_str(), markersize, ecolor.c_str(), markeredgecolor.c_str(), color.c_str());
  }
}

void matplotlibcpp::plot(vector<double> X, vector<double> Y, double lw, string color)
{
  if (X.size() == Y.size()) {
    fprintf(p, "plt.plot([");

    bool first = true;
    for (int i=0; i<X.size(); i++) {
      if (isfinite(X[i]) && isfinite(Y[i])) {
	if (first) {
	  fprintf(p,"%e",X[i]);
	  first = false;
	} else {
	  fprintf(p, ",%e", X[i]);
	}
      }
    }

    fprintf(p,"], [");

    first = true;
    for (int i=0; i<Y.size(); i++) {
      if (isfinite(X[i]) && isfinite(Y[i])) {
	if (first) {
	  fprintf(p,"%e",Y[i]);
	  first = false;
	} else {
	  fprintf(p, ",%e", Y[i]);
	}
      }
    }
    
    fprintf(p, "], lw=%f, color=\"%s\")\n", lw, color.c_str());
  }
}

void matplotlibcpp::xlog()
{
  fprintf(p, "plt.xscale(\"log\")\n");
}

void matplotlibcpp::ylog()
{
  fprintf(p, "plt.yscale(\"log\")\n");
}

void matplotlibcpp::contourf(vector<double> X, vector<double> Y, vector<double> Z, string label)
{
  int xsize = X.size();
  int ysize = Y.size();
  
  fprintf(p, "import numpy as np\n");
  fprintf(p, "x = np.array([%e", X[0]);
  for (int i=1; i<xsize; i++) {
    fprintf(p, ",%e", X[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "y = np.array([%e", Y[0]);
  for (int i=1; i<ysize; i++) {
    fprintf(p, ",%e", Y[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "X, Y = np.meshgrid(x,y)\n");
  fprintf(p, "Z = np.zeros((%d,%d))\n", ysize, xsize);
  for (int i=0; i<Z.size(); i++) {
    fprintf(p, "Z[%d,%d] = %e\n", i/xsize, i%xsize, Z[i]);
  }
  fprintf(p, "CF = plt.contourf(X,Y,Z)\n");
  fprintf(p, "CB = plt.colorbar(CF)\n");
  fprintf(p, "CB.set_label(\"%s\")\n", label.c_str());
}

void matplotlibcpp::log_contourf(vector<double> X, vector<double> Y, vector<double> Z,
				 string label)
{
  int xsize = X.size();
  int ysize = Y.size();
  
  fprintf(p, "import numpy as np\n");
  fprintf(p, "x = np.array([%e", X[0]);
  for (int i=1; i<xsize; i++) {
    fprintf(p, ",%e", X[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "y = np.array([%e", Y[0]);
  for (int i=1; i<ysize; i++) {
    fprintf(p, ",%e", Y[i]);
  }
  fprintf(p, "])\n");
  fprintf(p, "X, Y = np.meshgrid(x,y)\n");
  fprintf(p, "Z = np.zeros((%d,%d))\n", ysize, xsize);
  for (int i=0; i<Z.size(); i++) {
    fprintf(p, "Z[%d,%d] = %e\n", i/xsize, i%xsize, log10(Z[i]));
  }
  fprintf(p, "CF = plt.contourf(X,Y,Z)\n");
  fprintf(p, "CB = plt.colorbar(CF)\n");
  fprintf(p, "CB.set_label(\"%s\")\n", label.c_str());
}

void matplotlibcpp::show()
{
  fprintf(p, "plt.show()\n");
}

void matplotlibcpp::save(string filename)
{
  fprintf(p, "plt.savefig(\"%s\")\n", filename.c_str());
}

void matplotlibcpp::close()
{
  fprintf(p, "plt.close()");
  fprintf(p, "quit()");
}
  
#endif
