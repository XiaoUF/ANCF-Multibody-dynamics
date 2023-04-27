#include "Beam.h"
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_movstat.h>


class FSInterface
{
 public:
  // constructor
  FSInterface();
  FSInterface(Beam *, int);
  ~FSInterface();
  // allocate arrays used to transfer data between OpenFOAM and ANCF solver
  void setup(Beam *, int);
  void allocate(); 
  double **getPos();
  double **getPforce();
  double setPforce(double **, double **);
  // For particles
  void GetParticles(double *, double *, double *, double *, int);
  void GetBox(int, double &, double &, double &, double &);

  void calc_points();
  void calc_slopes();
  void Sme(double *, double *, double *);
  double distance(double, double, double);
  double distance(double *);
  void MaskU(double, double, double, double &, double *);
  double** getFs(); // interface to return fs to the structure solver
  double* getXs(); // interface to return xs to the structure solver
  double** getUc();
  double** getXl();
  double** getXr();
  double** getXc();
  double* getPl();
  double* getPr();
  double getWidth();
  double getThickness();
  void outputXUC();
  void hydroForce();
  void advanceRK2(double dt);
  void interpPressure(); // interpolate pressure from flow field
  void interpVelocity(); // interpolate velocity from flow field
  int getENP();
  int getTNP();
 private:
  Beam *bm;
  int ENP;
  int TNP;
  double **xc;
  double **sc;
  double **xl;
  double **xr;
  double *pl;
  double *pr;
  double *theta;
  double **pf; // force from the fluid solver
  double **uc; // velocity at the center
  // Point Force on the structure solver
  double **fs;
  double *xs;
  gsl_vector *ux;
  gsl_vector *uy;
  gsl_vector *uz;
  gsl_vector *uxm;
  gsl_vector *uym;
  gsl_vector *uzm;
  gsl_movstat_workspace *w;
};

