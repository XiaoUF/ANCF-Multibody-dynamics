#include "FSInterface.h"
#include <math.h>

FSInterface::FSInterface()
{
  bm  = NULL;
  ENP = 0;
  TNP = 0;
  xc = NULL;
  uc = NULL;
  xl = NULL;
  pr = NULL;
  pl = NULL;
  pf = NULL;
  theta = NULL;
  fs = NULL;
  xs = NULL;
  ux = NULL;
  uy = NULL;
  uz = NULL;
  uxm = NULL;
  uym = NULL;
  uzm = NULL;
  w   = NULL;
}

FSInterface::FSInterface(Beam *beam, int npe):
  bm(beam)
{
  ENP = npe;
  int ne = bm->getNE();
  TNP = npe+(ne-1)*(ENP-1);
  allocate();
}

void FSInterface::setup(Beam *beam, int npe)
{
  bm = beam;
  ENP = npe;
  int ne = bm->getNE();
  TNP = npe+(ne-1)*(ENP-1);
  allocate();
}

void FSInterface::allocate()
{
  int ne = bm->getNE();
  // Allocate memory
  int n1 = 3*TNP;
  double *data_xc, *data_xl, *data_xr;
  data_xc = (double *) malloc(sizeof(double)*n1);
  data_xl = (double *) malloc(sizeof(double)*n1);
  data_xr = (double *) malloc(sizeof(double)*n1);
  int n=0;
  xc = (double **) malloc(3*sizeof(double *));
  xr = (double **) malloc(3*sizeof(double *));
  xl = (double **) malloc(3*sizeof(double *));
  for(int i=0; i<3; i++){
    xc[i] = &data_xc[n];
    xr[i] = &data_xr[n];
    xl[i] = &data_xl[n];
    n = n+TNP;
  }

  double *data_sc;
  data_sc = (double *) malloc(sizeof(double)*n1);
  sc = (double **) malloc(3*sizeof(double *));
  n = 0;
  for(int i=0; i<3; i++){
    sc[i] = &data_sc[n];
    n = n+TNP;
  }

  double *data_pf, *data_uc;
  data_pf = (double *) malloc(sizeof(double)*n1);
  data_uc = (double *) malloc(sizeof(double)*n1);
  pf = (double **) malloc(3*sizeof(double *));
  uc = (double **) malloc(3*sizeof(double *));
  n = 0;
  for(int i=0; i<3; i++){
    pf[i] = &data_pf[n];
    uc[i] = &data_uc[n];
    n = n+TNP;
  }

  pl = (double *) malloc(sizeof(double)*TNP);
  pr = (double *) malloc(sizeof(double)*TNP);
  theta = (double *) malloc(sizeof(double)*TNP);

  // Allocate the data to the structure solver
  int nn = 3*ne;
  xs = (double *) malloc(sizeof(double)*ne);
  double *data_fs;
  data_fs = (double *) malloc(sizeof(double)*nn);
  n = 0;
  fs = (double **) malloc(sizeof(double *)*3);
  for(int i=0; i<3; i++){
    fs[i] = &data_fs[n];
    n = n+ne;
  }

  // Allocate data
  const size_t K = 5;
  ux  = gsl_vector_alloc(ne);
  uxm = gsl_vector_alloc(ne);
  uy  = gsl_vector_alloc(ne);
  uym = gsl_vector_alloc(ne);
  uz  = gsl_vector_alloc(ne);
  uzm = gsl_vector_alloc(ne);
  w = gsl_movstat_alloc(K);
}

FSInterface::~FSInterface()
{
  free(pl);
  free(pr);
  free(theta);
  free(xc[0]);
  free(xl[0]);
  free(xr[0]);
  free(sc[0]);
  free(pf[0]);
  free(uc[0]);
  free(fs[0]);
  free(xc);
  free(xl);
  free(xr);
  free(pf);
  free(fs);
  free(xs);
  free(uc);
  gsl_vector_free(ux);
  gsl_vector_free(uy);
  gsl_vector_free(uz);
  gsl_vector_free(uxm);
  gsl_vector_free(uym);
  gsl_vector_free(uzm);
  gsl_movstat_free(w);
}

void FSInterface::GetParticles(double *xp, double *zp, double *up, double *wp, int npe)
{ 
  int ie;
  int ip;
  int ii = 0;
  double ee[12];
  double uu[12];
  double S[4];
  double x,eta;
  double elen;
  gsl_vector *eb;
  gsl_vector *eu;
  int ne;
  
  ne = bm->getNE();
  eb = bm->getEE();
  eu = bm->getU();
  elen = bm->getLength()/double(ne);
  
  for(ie=0; ie<ne; ie++){
    int istart = 6*ie; 
    for(int j=0; j<12; j++){
      ee[j] = gsl_vector_get(eb,istart+j);
      uu[j] = gsl_vector_get(eu,istart+j);
    }
    for(int ip=0; ip<npe; ip++){
      eta = (double(ip)+0.5)/double(npe);
      S[0] = 1.-3.*eta*eta + 2.0*eta*eta*eta;
      S[1] = elen*(eta-2.*eta*eta+eta*eta*eta);
      S[2] = 3.*eta*eta-2.*eta*eta*eta;
      S[3] = elen*(-eta*eta+eta*eta*eta);
      xp[ii] = S[0]*ee[0] + S[1]*ee[3] + S[2]*ee[6] + S[3]*ee[9];
      zp[ii] = S[0]*ee[2] + S[1]*ee[5] + S[2]*ee[8] + S[3]*ee[11];
      up[ii] = S[0]*uu[0] + S[1]*uu[3] + S[2]*uu[6] + S[3]*uu[9];
      wp[ii] = S[0]*uu[2] + S[1]*uu[5] + S[2]*uu[8] + S[3]*uu[11];
      ii++;
    }
  }
}

void FSInterface::GetBox(int ne, double &xlo, double &zlo, double &xhi, double &zhi)
{ 
  gsl_vector *eb;
  eb = bm->getEE();
  double x1, x2;
  double z1, z2;
  int istart;
  istart = 6*ne;
  x1 = gsl_vector_get(eb,istart);
  x2 = gsl_vector_get(eb,istart+6);
  z1 = gsl_vector_get(eb,istart+2);
  z2 = gsl_vector_get(eb,istart+8);
  if(x1>x2){
    xhi = x1;
    xlo = x2;
  } else{ 
    xhi = x2;
    xlo = x1;
  }
  
  if(z1>z2){
    zlo = z2;
    zhi = z1;
  } else{ 
    zhi = z2;
    zlo = z1;
  }

}

void FSInterface::calc_points()
{
  int ie; 
  int ip;
  int ii = 0;
  double ee[12];
  double uu[12];
  double S[4];
  double x,eta;
  double elen;
  gsl_vector *eb;
  gsl_vector *eu;
  int ne;

  ne = bm->getNE();
  eb = bm->getEE();
  eu = bm->getU();
  elen = bm->getLength()/double(ne);

// filter the velocity
  for(ie=0; ie<ne; ie++){
     gsl_vector_set(ux,ie,gsl_vector_get(eu,6*ie));
     gsl_vector_set(uy,ie,gsl_vector_get(eu,6*ie+1));
     gsl_vector_set(uz,ie,gsl_vector_get(eu,6*ie+2));
  }
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,ux,uxm,w);
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,uy,uym,w);
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,uz,uzm,w);

  for(ie=1; ie<ne; ie++){
     gsl_vector_set(eu,6*ie,gsl_vector_get(uxm,ie));
     gsl_vector_set(eu,6*ie+1,gsl_vector_get(uym,ie));
     gsl_vector_set(eu,6*ie+2,gsl_vector_get(uzm,ie));
  }
// Position
/*
  for(ie=0; ie<ne; ie++){
     gsl_vector_set(ux,ie,gsl_vector_get(eb,6*ie));
     gsl_vector_set(uy,ie,gsl_vector_get(eb,6*ie+1));
     gsl_vector_set(uz,ie,gsl_vector_get(eb,6*ie+2));
  }
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,ux,uxm,w);
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,uy,uym,w);
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,uz,uzm,w);

  for(ie=1; ie<ne; ie++){
     gsl_vector_set(eb,6*ie,gsl_vector_get(uxm,ie));
     gsl_vector_set(eb,6*ie+1,gsl_vector_get(uym,ie));
     gsl_vector_set(eb,6*ie+2,gsl_vector_get(uzm,ie));
  }
// The orientation
// */
  for(ie=0; ie<ne; ie++){
    int istart = 6*ie;
    // get the length of the element
    for(int j=0; j<12; j++){
      ee[j] = gsl_vector_get(eb,istart+j);
      uu[j] = gsl_vector_get(eu,istart+j);
    }
    for(int ip=0; ip<ENP-1; ip++){
      // get the shape function
      eta = double(ip)/double(ENP);
      S[0] = 1.-3.*eta*eta + 2.0*eta*eta*eta;
      S[1] = elen*(eta-2.*eta*eta+eta*eta*eta);
      S[2] = 3.*eta*eta-2.*eta*eta*eta;
      S[3] = elen*(-eta*eta+eta*eta*eta);
      xc[0][ii] = S[0]*ee[0] + S[1]*ee[3] + S[2]*ee[6] + S[3]*ee[9];
      xc[1][ii] = S[0]*ee[1] + S[1]*ee[4] + S[2]*ee[7] + S[3]*ee[10];
      xc[2][ii] = S[0]*ee[2] + S[1]*ee[5] + S[2]*ee[8] + S[3]*ee[11];
      uc[0][ii] = S[0]*uu[0] + S[1]*uu[3] + S[2]*uu[6] + S[3]*uu[9];
      uc[1][ii] = S[0]*uu[1] + S[1]*uu[4] + S[2]*uu[7] + S[3]*uu[10];
      uc[2][ii] = S[0]*uu[2] + S[1]*uu[5] + S[2]*uu[8] + S[3]*uu[11];
      ii++;
    }
  }
  // End point
  xc[0][ii] = ee[6];
  xc[1][ii] = ee[7];
  xc[2][ii] = ee[8];
  uc[0][ii] = uu[6];
  uc[1][ii] = uu[7];
  uc[2][ii] = uu[8];
}

// Debug use only
void FSInterface::outputXUC()
{
  int i;
  char fname[32] = "xuc.txt";
  FILE *fp;
  fp = fopen(fname,"w+");
  for(i=0; i<TNP; i++){
    fprintf(fp,"%12.9e\t%12.9e\t%12.9e\t%12.9e\t%12.9e\t%12.9e\n",xc[0][i],xc[1][i],xc[2][i],uc[0][i],uc[1][i],uc[2][i]);
  }
  fclose(fp);
}


void FSInterface::calc_slopes()
{
  int ie; 
  int ip;
  int ii = 0;
  double ee[12];
  double S[4];
  double x,eta;
  double elen;
  gsl_vector *eb;
  int ne;

  ne = bm->getNE();
  eb = bm->getEE();
  elen = bm->getLength()/double(ne);

  for(ie=0; ie<ne; ie++){
    int istart = 6*ie;
    // get the length of the element
    for(int j=0; j<12; j++)
      ee[j] = gsl_vector_get(eb,istart+j);
    for(ip=0; ip<ENP-1; ip++){
      // get the shape function
      eta = double(ip)/double(ENP);
      S[0] = (6.*eta*eta - 6*eta)/elen;
      S[1] = 1.-4.*eta+3.*eta*eta;
      S[2] = (-6.*eta*eta+6.*eta)/elen;
      S[3] = -2.*eta+3.*eta*eta;
      sc[0][ii] = S[0]*ee[0] + S[1]*ee[3] + S[2]*ee[6] + S[3]*ee[9];;
      sc[1][ii] = S[0]*ee[1] + S[1]*ee[4] + S[2]*ee[7] + S[3]*ee[10];
      sc[2][ii] = S[0]*ee[2] + S[1]*ee[5] + S[2]*ee[8] + S[3]*ee[11];
      ii++;
    }
  }
  // End point
  sc[0][ii] = ee[9];
  sc[1][ii] = ee[10];
  sc[2][ii] = ee[11];

  for(int i=0; i<TNP; i++)
    theta[i] = atan2(sc[2][i],sc[0][i]);

  double thickness = bm->getThickness()*0.5;
  for(int i=0; i<TNP; i++){
    xr[0][i] = xc[0][i] + thickness*sin(theta[i]);
    xr[1][i] = xc[1][i];
    xr[2][i] = xc[2][i] - thickness*cos(theta[i]);
    xl[0][i] = xc[0][i] - thickness*sin(theta[i]);
    xl[1][i] = xc[1][i];
    xl[2][i] = xc[2][i] + thickness*cos(theta[i]);
  }

}

void FSInterface::MaskU(double x, double y, double z, double &dist, double *um)
{ 
  int i, ind; 
  double d2;
  double dx,dy,dz;
  dist = 1.e10;
  ind  = 0;
  for(i=0; i<TNP; i++){
    dx = x-xc[0][i];
    dz = z-xc[2][i];
    d2 = sqrt(dx*dx+dz*dz);
    if(d2<dist) {
       dist = d2;
       ind = i;
    }
  }
  um[0] = uc[0][ind];
  um[1] = uc[1][ind];
  um[2] = uc[2][ind];
}

double FSInterface::distance(double x, double y, double z)
{
  int i;
  double dist;
  double d2;
  double dx,dy,dz;
  dist = 1.e10;
  for(i=0; i<TNP; i++){
    dx = x-xc[0][i];
    //dy = y-xc[1][i];
    dz = z-xc[2][i];
    d2 = sqrt(dx*dx+dz*dz);
    if(d2<dist) 
      dist = d2;
  }
  return dist;
}

double FSInterface::distance(double *xp)
{
  int i;
  double dist;
  double d2;
  double dx,dy,dz;
  dist = 1.e10;
  for(i=0; i<TNP; i++){
    dx = xp[0]-xc[0][i];
    dz = xp[2]-xc[2][i];
    d2 = sqrt(dx*dx+dz*dz);
    if(d2<dist) 
      dist = d2;
  }
  return dist;

}

void FSInterface::hydroForce()
{
  int i;
  double len = bm->getLength();

//  double dpdz = 200.;
#if 0
  for(i=0; i<TNP; i++){
    pl[i] = dpdz*xl[2][i];;
    pr[i] = dpdz*xr[2][i];
  }
#endif
  
  for(i=0; i<TNP; i++){
    pf[0][i] = (pl[i]-pr[i])*sin(theta[i]);
    pf[1][i] = 0.0;
    pf[2][i] = -(pl[i]-pr[i])*cos(theta[i]);
  }  

  int ne = bm->getNE();
  double alpha = len/double(ne)/double(ENP-1);
  double beta  = alpha*bm->getWidth();
#if 0
  std::cout << "length = " << bm->getLength() << std::endl;
  std::cout << "width = " << bm->getWidth() << std::endl;
  std::cout << "alpha = " << alpha << std::endl;
  std::cout << "beta = " << beta << std::endl;
#endif
  // get the center of force
  double spe;
  double sp;
  double dp1,dp2,dp;
  for(i=0; i<ne; i++){
    spe = 0.;
    sp  = 0.;
    int istart = (ENP-1)*i;
    for(int j=1; j<ENP; j++){
      dp2 = pl[istart+j] - pr[istart+j];
      dp1 = pl[istart+j-1] - pr[istart+j-1];
      dp = 0.5*(dp1+dp2);
      sp = sp+dp;
      spe = spe+dp*(double(j)-0.5);
    }
    //std::cout << "XS = " << alpha*spe/sp << std::endl;
    if(fabs(sp) > 1.e-8){
      xs[i] = alpha*spe/sp;
     // std::cout << "xs = " << xs[i] << std::endl;
    }
    else
      xs[i] = 0.;
  }
  // sum over the force
  for(i=0; i<ne; i++){
    int istart;
    int iend;
    istart = (ENP-1)*i;
    iend   = istart + ENP-1;
    fs[0][i] = 0.0;
    fs[1][i] = 0.0;
    fs[2][i] = 0.0;
    for(int j=istart+1; j<iend; j++){
      fs[0][i] = fs[0][i] + pf[0][j];
      fs[1][i] = fs[1][i] + pf[1][j];
      fs[2][i] = fs[2][i] + pf[2][j];
    }  
    fs[0][i] = fs[0][i] + 0.5*(pf[0][istart]+pf[0][iend]);
    fs[1][i] = fs[1][i] + 0.5*(pf[1][istart]+pf[1][iend]);
    fs[2][i] = fs[2][i] + 0.5*(pf[2][istart]+pf[2][iend]);

    // scale
    fs[0][i] = fs[0][i]*beta;
    fs[1][i] = fs[1][i]*beta;
    fs[2][i] = fs[2][i]*beta;
  }

  // Smooth the force
#if 0
  for(i=0; i<ne; i++){
     gsl_vector_set(ux,i,fs[0][i]);
     gsl_vector_set(uy,i,fs[1][i]);
     gsl_vector_set(uz,i,fs[2][i]);
  }
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,ux,uxm,w);
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,uy,uym,w);
  gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,uz,uzm,w);
  for(i=0; i<ne; i++){
     fs[0][i] = gsl_vector_get(uxm,i);
     fs[1][i] = gsl_vector_get(uym,i);
     fs[2][i] = gsl_vector_get(uzm,i);
  }
#endif
}

double ** FSInterface::getFs()
{
  return fs;
}

double * FSInterface::getXs()
{
  return xs;
}

double ** FSInterface::getXc()
{
  return xc;
}

double ** FSInterface::getUc()
{
  return uc;
}

double ** FSInterface::getXl()
{
  return xl;
}

double ** FSInterface::getXr()
{
  return xr;
}

double * FSInterface::getPl()
{
  return pl;
}

double * FSInterface::getPr()
{
  return pr;
}

int FSInterface::getTNP()
{
  return TNP;
}

double FSInterface::getWidth()
{
  return bm->getWidth();
}

double FSInterface::getThickness()
{ 
  return bm->getThickness();
}
