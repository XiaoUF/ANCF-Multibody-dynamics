#ifndef BEAM_H
#define BEAM_H

#include "node.h"
#include "Element.h"
#include "particle.h"
#include "Gauss.h"
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#define DEBUG 1

class Beam{
 private:
  int ne;
  double length;
  double width;
  double thickness;
  double E;
  double nu;
  double rho;
  double mass;
  // For the elements
  Element *eb;
  // Gauss points
  GaussPoint gs;
  // Position, velocity and force vectors
  double **v;      
  double **x;
  double **f;
  // Mass matrix
  gsl_matrix *MM;
  gsl_matrix *MB;
  gsl_matrix *Me;
  gsl_vector_uint *piv;
  gsl_vector *rhs;
  gsl_vector *sol;
  gsl_vector *ee; 
  gsl_vector *eo;
  gsl_vector *et; // temporary vector holder
  gsl_vector *u;  // Velocity vector
  gsl_vector *uo;
  gsl_vector *ut;
  // Force vectors
  gsl_vector *qg;
  gsl_vector *qh;
  gsl_vector *qe;
  gsl_vector *qc;

  
 public:
  // Default constructor
  Beam(){
    ne = 0;
    length = 0.0;
    width = 0.0;
    thickness = 0.0;
    E = 0.0;
    nu = 0.0;
    rho = 0.0;
    mass = 0.0;
    // Pointers
    eb = NULL;
    v = NULL;
    x = NULL;
    f = NULL;
    MM = NULL;
    MB = NULL;
    Me = NULL;
    piv = NULL;
    rhs = NULL;
    sol = NULL;
    ee = NULL;
    eo = NULL;
    et = NULL;
    u = NULL;
    uo = NULL;
    ut = NULL;
    qg = NULL;
    qh = NULL;
    qe = NULL;
    qc = NULL;
  }
  
  // Constructor
  Beam(int ne, double length, vec3 x0, int dir, double width, double thickness){
    // set-up the nodes
    this->ne = ne;

    // Memory management
    allocate();
    this->length = length;
    this->width  = width;
    this->thickness = thickness;
    double dl = length/ne;
    int nd = ne+1;
    // Direction vector
    double dv[3] = {0., 0., 0.};
    dv[dir] = 1.0;
    
    switch(dir){
    case 0:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x + i*dl;
	x[1][i] = x0.y;
	x[2][i] = x0.z;
      }
      break;
    case 1:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x;
	x[1][i] = x0.y + i*dl;
	x[2][i] = x0.z;
      }
      break;
    case 2:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x;
	x[1][i] = x0.y;
	x[2][i] = x0.z + i*dl;
      }
      break;
    default:
      std::cout << "Dir should be in the range of [0,2]" << std::endl;
    }

    // Setup the Elements
    for(int i=0; i<ne; i++){
      Node node0 = Node(x[0][i],x[1][i],x[2][i],dv[0],dv[1],dv[2]);
      Node node1 = Node(x[0][i+1],x[1][i+1],x[2][i+1],dv[0],dv[1],dv[2]);
      eb[i] = Element(node0, node1, width, thickness);
    }
    AssembleMassMatrix();
    InitialCondition(x0,dir);
  }

  Beam(int ne, double length, vec3 x0, int dir, double width, double thickness, double E, double nu, double rho){
    // set-up the nodes
    this->ne = ne;
    allocate();
    this->length = length;
    this->width  = width;
    this->thickness = thickness;
    this->E   = E;
    this->nu  = nu;
    this->rho = rho;
    
    double dl = length/ne;
    int nd = ne+1;
    // Direction vector
    double dv[3] = {0., 0., 0.};
    dv[dir] = 1.0;
    
    switch(dir){
    case 0:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x + i*dl;
	x[1][i] = x0.y;
	x[2][i] = x0.z;
      }
      break;
    case 1:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x;
	x[1][i] = x0.y + i*dl;
	x[2][i] = x0.z;
      }
      break;
    case 2:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x;
	x[1][i] = x0.y;
	x[2][i] = x0.z + i*dl;
      }
      break;
    default:
      std::cout << "Dir should be in the range of [0,2]" << std::endl;
    }
    // Setup the Elements
    for(int i=0; i<ne; i++){
      Node node0 = Node(x[0][i],x[1][i],x[2][i],dv[0],dv[1],dv[2]);
      Node node1 = Node(x[0][i+1],x[1][i+1],x[2][i+1],dv[0],dv[1],dv[2]);
      eb[i] = Element(node0, node1, width, thickness, E, nu, rho);
    }
    AssembleMassMatrix();
    InitialCondition(x0, dir);
  }

  // Destructor
  ~Beam() {
    // Deallocate memory
    free(eb);
    // Deallocate arrays
    free(v[0]);
    free(x[0]);
    free(f[0]);
    free(x);
    free(f);
    free(v);

    // Deallocate the mass matrix
    gsl_matrix_free(MM);
    gsl_matrix_free(MB);
    gsl_matrix_free(Me);
    gsl_vector_free(rhs);
    gsl_vector_free(sol);
    gsl_vector_free(ee);
    gsl_vector_free(eo);
    gsl_vector_free(et);
    gsl_vector_free(u);
    gsl_vector_free(uo);
    gsl_vector_free(ut);
    gsl_vector_free(qg);
    gsl_vector_free(qh);
    gsl_vector_free(qe);
    gsl_vector_free(qc);
    gsl_vector_uint_free(piv);
  }

  void allocate(){
    // use 2D array
    int nd,n1,n2;
    nd = ne+1;
    n1 = nd*3;
    n2 = 6*nd;
    // Allocate data for the data
    double *data_v, *data_p, *data_f;
    data_p = (double *) malloc(sizeof(double)*n1);
    data_v = (double *) malloc(sizeof(double)*n1);
    data_f = (double *) malloc(sizeof(double)*n1);
    v = (double **) malloc(3*sizeof(double *));
    x = (double **) malloc(3*sizeof(double *));
    f = (double **) malloc(3*sizeof(double *));
    int n = 0;
    for(int i=0; i<3; i++){
      v[i] = &data_v[n];
      x[i] = &data_p[n];
      f[i] = &data_f[n];
      n = n+nd;
    }

    // Initialize the arrays
    for(int i=0; i<nd; i++){
      v[0][i] = 0.0;
      v[1][i] = 0.0;
      v[2][i] = 0.0;
      x[0][i] = 0.0;
      x[1][i] = 0.0;
      x[2][i] = 0.0;
      f[0][i] = 0.0;
      f[1][i] = 0.0;
      f[2][i] = 0.0;
    }

    // Allocate memory for elements in the beam
    eb = (Element *) malloc(ne*sizeof(Element));

    // Allocate data for the mass matrix
    int p = 11;
    int q = 11;
    int nb = 2*q+p+1;
    MM = gsl_matrix_alloc(n2,n2);
    MB = gsl_matrix_alloc(n2,nb);
    Me = gsl_matrix_alloc(12,12);
    piv = gsl_vector_uint_alloc(n2);
    qg  = gsl_vector_alloc(n2);
    qe  = gsl_vector_alloc(n2);
    qh  = gsl_vector_alloc(n2);
    qc  = gsl_vector_alloc(n2);
    rhs = gsl_vector_alloc(n2);
    sol = gsl_vector_alloc(n2);
    ee = gsl_vector_alloc(n2);
    eo = gsl_vector_alloc(n2);
    et = gsl_vector_alloc(n2);
    u  = gsl_vector_alloc(n2);
    uo = gsl_vector_alloc(n2);
    ut = gsl_vector_alloc(n2);
    
    // Initialize the matrix
    gsl_matrix_set_zero(MM);
    gsl_matrix_set_zero(MB);
    gsl_vector_set_zero(sol);
    gsl_vector_set_zero(rhs);
    gsl_vector_set_zero(qg);
    gsl_vector_set_zero(qe);
    gsl_vector_set_zero(qh);
    gsl_vector_set_zero(qc);
    gsl_vector_set_zero(ee);
    gsl_vector_set_zero(eo);
    gsl_vector_set_zero(et);
    gsl_vector_set_zero(u);
    gsl_vector_set_zero(uo);
  }

  void SetupBeam(int ne, double length, vec3 x0, int dir, double width, double thickness, double E, double nu, double rho){
   // set-up the nodes
    this->ne = ne;
    allocate();
    this->length = length;
    this->width  = width;
    this->thickness = thickness;
    this->E   = E;
    this->nu  = nu;
    this->rho = rho;
    
    double dl = length/ne;
    int nd = ne+1;
    // Direction vector
    double dv[3] = {0., 0., 0.};
    dv[dir] = 1.0;
    
    switch(dir){
    case 0:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x + i*dl;
	x[1][i] = x0.y;
	x[2][i] = x0.z;
      }
      break;
    case 1:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x;
	x[1][i] = x0.y + i*dl;
	x[2][i] = x0.z;
      }
      break;
    case 2:
      for(int i=0; i<nd; i++){
	x[0][i] = x0.x;
	x[1][i] = x0.y;
	x[2][i] = x0.z + i*dl;
      }
      break;
    default:
      std::cout << "Dir should be in the range of [0,2]" << std::endl;
    }
    // Setup the Elements
    for(int i=0; i<ne; i++){
      Node node0 = Node(x[0][i],x[1][i],x[2][i],dv[0],dv[1],dv[2]);
      Node node1 = Node(x[0][i+1],x[1][i+1],x[2][i+1],dv[0],dv[1],dv[2]);
      eb[i] = Element(node0, node1, width, thickness, E, nu, rho);
    }
    AssembleMassMatrix();
    InitialCondition(x0, dir);
  }

  int gen2band_matrix(const size_t p, const size_t q, const gsl_matrix * A, gsl_matrix * AB)
  {
    const size_t N = A->size2;

    if (AB->size1 != N)
      {
	GSL_ERROR("banded matrix requires N rows", GSL_EBADLEN);
      }
    else if (AB->size2 != 2*p + q + 1)
      {
	GSL_ERROR("banded matrix requires 2*p + q + 1 columns", GSL_EBADLEN);
      }
    else
      {
	size_t i;

	gsl_matrix_set_zero(AB);

	/* copy diagonal and subdiagonals */
	for (i = 0; i <= p; ++i)
	  {
	    gsl_vector_const_view v = gsl_matrix_const_subdiagonal(A, i);
	    gsl_vector_view w = gsl_matrix_subcolumn(AB, p + q + i, 0, v.vector.size);
	    gsl_vector_memcpy(&w.vector, &v.vector);
	  }

	/* copy superdiagonals */
	for (i = 1; i <= q; ++i)
	  {
	    gsl_vector_const_view v = gsl_matrix_const_superdiagonal(A, i);
	    gsl_vector_view w = gsl_matrix_subcolumn(AB, p + q - i, i, v.vector.size);
	    gsl_vector_memcpy(&w.vector, &v.vector);
	  }

	return GSL_SUCCESS;
      }
  }

  // setup mass matrix
  void AssembleMassMatrix()
  {
    int s;
    int p = 11;
    int q = 11;
    int n1 = 6*(ne+1);
    // Setup the mass matrix
    //std::cout << "Setup the mass matrix" << std::endl;
    for(size_t i=0; i<ne; i++){
      ComputeMassMatrix(eb[i], Me);
      double coef = eb[i].getDensity()*eb[i].getArea()*eb[i].getLength_l();
      gsl_matrix_scale(Me,coef);
      //std::cout << "Me[0][0] = " << gsl_matrix_get(Me,0,0) << std::endl;
      size_t ii = i*6;
      for(size_t j=0; j<6; j++)
	for(size_t k=0; k<6; k++){
	  gsl_matrix_set(MM,ii+j,ii+k,gsl_matrix_get(MM,ii+j,ii+k)+gsl_matrix_get(Me,j,k));
	  gsl_matrix_set(MM,ii+j,ii+k+6,gsl_matrix_get(Me,j,k+6));
	  gsl_matrix_set(MM,ii+j+6,ii+k,gsl_matrix_get(Me,j+6,k));
	  gsl_matrix_set(MM,ii+j+6,ii+k+6,gsl_matrix_get(MM,ii+j+6,ii+k+6)+gsl_matrix_get(Me,j+6,k+6));
	}
    }

    // Setup the constraint at the root
    for(size_t j=0; j<6; j++)
      for(size_t k=0; k<n1; k++){
	gsl_matrix_set(MM,j,k,0.0);
        gsl_matrix_set(MM,k,j,0.0);
    }
    for(size_t j=0; j<6; j++)
      gsl_matrix_set(MM,j,j,1.0);

    // construct band matrix
    gen2band_matrix(p,q,MM,MB);

    // LU decomposition with pivoting
    s = gsl_linalg_LU_band_decomp(n1,p,q,MB,piv);

  }

  void InitialCondition(vec3 x0, int dir)
  {
    // Initialize ee and eo vectors
    int i;
    int nd = ne+1;
    double dl = length/ne;
    switch(dir){
    case 0:
      for(i=0; i<nd; i++){
	gsl_vector_set(ee,6*i,x0.x+i*dl);
	gsl_vector_set(ee,6*i+1,x0.y);
	gsl_vector_set(ee,6*i+2,x0.z);
	gsl_vector_set(ee,6*i+3,1.0);
      }
      break;
    case 1:
      for(i=0; i<nd; i++){
	gsl_vector_set(ee,6*i,x0.x);
	gsl_vector_set(ee,6*i+1,x0.y+i*dl);
	gsl_vector_set(ee,6*i+2,x0.z);
	gsl_vector_set(ee,6*i+4,1.0);
      }
      break;
    case 2:
      for(i=0; i<nd; i++){
	gsl_vector_set(ee,6*i,x0.x);
	gsl_vector_set(ee,6*i+1,x0.y);
	gsl_vector_set(ee,6*i+2,x0.z+i*dl);
	gsl_vector_set(ee,6*i+5,1.0);
      }
      break;
    default:
      std::cout << "Dir should be in the range of [0,2]" << std::endl;
    }

#if 0
    FILE *fp;
    fp=fopen("ee.txt","w+");
    gsl_vector_fprintf(fp,ee,"%9.6e");
    fclose(fp);
#endif
  }

  void ComputeMassMatrix(Element el, gsl_matrix *Me)
  {
    int i;
    double l = el.getLength_l();
    double t2 = l*11./210.;
    double t3 = l*l;
    double t4 = t3/105.;
    double t5 = l*13./420.;

    gsl_matrix_set_zero(Me);

    for(i=0;i<3;i++){
      // First 3 rows
      gsl_matrix_set(Me,i,i,13./35.);
      gsl_matrix_set(Me,i,i+3,t2);
      gsl_matrix_set(Me,i,i+6,9./70.);
      gsl_matrix_set(Me,i,i+9,-t5);
      // Second 3 rows
      gsl_matrix_set(Me,3+i,i,t2);
      gsl_matrix_set(Me,3+i,3+i,t4);
      gsl_matrix_set(Me,3+i,6+i,t5);
      gsl_matrix_set(Me,i+3,i+9,-t3/140.);
      // Third 3 rows
      gsl_matrix_set(Me,6+i,i,9./70.);
      gsl_matrix_set(Me,6+i,i+3,t5);
      gsl_matrix_set(Me,i+6,i+6,13./35.);
      gsl_matrix_set(Me,6+i,i+9,-t2);
      //
      gsl_matrix_set(Me,9+i,i,-13.*l/420.);
      gsl_matrix_set(Me,i+9,i+3,-t3/140.);
      gsl_matrix_set(Me,i+9,i+6,-t2);
      gsl_matrix_set(Me,i+9,i+9,t4);
    }

  }

  void eps_integrand(int ie, double *e, double x, double *res)
  {
    double *Sx;
    Sx = eb[ie].shapeFunc(x,1);
    eb[ie].eps_eps_e(Sx, e, res);
  }

  void kappa_integrand(int ie, double *e, double x, double *res)
  {
    double Sx[4];
    double Sxx[4];
    double rx[3];
    double rxx[3];
    double v[3];
    double coef;
    double inv_rx_rx_cubed;
    eb[ie].shapeFunc(x,1,&Sx[0]);
    eb[ie].shapeFunc(x,2,&Sxx[0]);
    eb[ie].S_e(&Sx[0], e, &rx[0]);
    eb[ie].S_e(&Sxx[0], e, &rxx[0]);
    v[0] = rx[1]*rxx[2]-rx[2]*rxx[1];
    v[1] = rx[2]*rxx[0]-rx[0]*rxx[2];
    v[2] = rx[0]*rxx[1]-rx[1]*rxx[0];
    double rx3;
    rx3 = rx[0]*rx[0]+rx[1]*rx[1]+rx[2]*rx[2];
    coef = 3.*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/rx3;
    inv_rx_rx_cubed = 1.0/(rx3*rx3*rx3);
    eb[ie].kappa_kappa_e(&Sx[0], &Sxx[0], coef, inv_rx_rx_cubed, &rx[0], &rxx[0], &v[0], res);
  }

  void elasticForces()
  {
    int i,j,k;
    gsl_vector_set_zero(qe);
    for(i=0; i<ne; i++){
      double Qe[12];
      double res[12];
      double e[12];
      double coefe = 0.5*eb[i].getYoungsModulus()*eb[i].getLength_l()*eb[i].getArea();
      double coefk = 0.5*eb[i].getYoungsModulus()*eb[i].getLength_l()*eb[i].getMomentInertia();
      // Copy data from vector ee
      for(j=0; j<12; j++)
	e[j] = gsl_vector_get(ee,6*i+j);
      // Calculate axial force
      for (j=0; j<12; j++)
	Qe[j] = 0.0;
      for(j=0; j<5; j++){
	double x = (1.0+gs.x5[j])*eb[i].getLength_l()*0.5;
	eps_integrand(i,&e[0],x,&res[0]);
	for(k=0; k<12; k++){
	  Qe[k] = Qe[k] + coefe*gs.w5[j]*res[k];
	}
      }
      // calculate torsional force
      double Qk[12];
      for(j=0; j<12; j++)
	Qk[j] = 0.0;
      for(j=0; j<3; j++){
	double x = (1.0+gs.x3[j])*eb[i].getLength_l()*0.5;
	kappa_integrand(i, &e[0], x, res);
	for(k=0; k<12; k++){
	  Qk[k] = Qk[k] + coefk*gs.w3[j]*res[k];
	}
      }
      // Setup the elastic force
      for(j=0; j<12; j++)
	gsl_vector_set(qe,6*i+j, gsl_vector_get(qe,6*i+j)+Qe[j]+Qk[j]); 
    }
     // Boundary condition - fix the first node
     for(i=0; i<6; i++)
        gsl_vector_set(qe,i,0.0);
  }

  void pointLoad(int ie, double x, double *F)
  {
    double *Sx;
    Sx = eb[ie].shapeFunc(x,0);
    double Qp[12];
    Qp[0]  = F[0]*Sx[0];
    Qp[1]  = F[1]*Sx[0];
    Qp[2]  = F[2]*Sx[0];
    Qp[3]  = F[0]*Sx[1];
    Qp[4]  = F[1]*Sx[1];
    Qp[5]  = F[2]*Sx[1];
    Qp[6]  = F[0]*Sx[2];
    Qp[7]  = F[1]*Sx[2];
    Qp[8]  = F[2]*Sx[2];
    Qp[9]  = F[0]*Sx[3];
    Qp[10] = F[1]*Sx[3];
    Qp[11] = F[2]*Sx[3];

    int istart = 6*ie;
    for(int i=0; i<12; i++)
      gsl_vector_set(qh,istart+i,gsl_vector_get(qh,istart+i)+Qp[i]);
  }
  
  void gravityForce(double *g)
  {
    int i;
    gsl_vector_set_zero(qg);
    for(i=0; i<ne; i++){
      double l = eb[i].getLength_l();
      double t2 = g[0];
      double t3 = g[1];
      double t4 = g[2];
      double t5 = 0.5*t2;
      double t6 = 0.5*t3;
      double t7 = 0.5*t4;
      double t8 = l*t2/12.;
      double t9 = l*t3/12.;
      double t10 = l*t4/12.;
      double Qg[12] = {t5, t6, t7, t8, t9, t10, t5, t6, t7, -t8, -t9, -t10};
      double coef = eb[i].getDensity()*l*eb[i].getWidth()*eb[i].getThickness();
      for(int j=0; j<12; j++)
	gsl_vector_set(qg,6*i+j,gsl_vector_get(qg,6*i+j)+Qg[j]*coef);
    }
  }

  void setBeam(double *data)
  {
    
    for(int i=0; i<ne; i++){
      int i0 = i*6;
      int i1 = (i+1)*6;
      Node node0 = Node(data[i0],data[i0+1],data[i0+2],data[i0+3],data[i0+4],data[i0+5]);
      Node node1 = Node(data[i1],data[i1+1],data[i1+2],data[i1+3],data[i1+4],data[i1+5]);
      eb[i].UpdateNodes(node0,node1);
      //eb[i].printInfo();
    }
  }

  void updateBeam()
  {
    // update element using values from the updated vector ee
    for(int i=0; i<ne; i++){
      int i0 = i*6;
      int i1 = (i+1)*6;
      double x,y,z,dx,dy,dz;
      x = gsl_vector_get(ee,i0);
      y = gsl_vector_get(ee,i0+1);
      z = gsl_vector_get(ee,i0+2);
      dx = gsl_vector_get(ee,i0+3);
      dy = gsl_vector_get(ee,i0+4);
      dz = gsl_vector_get(ee,i0+5);
      Node node0 = Node(x,y,z,dx,dy,dz);
      x = gsl_vector_get(ee,i1);
      y = gsl_vector_get(ee,i1+1);
      z = gsl_vector_get(ee,i1+2);
      dx = gsl_vector_get(ee,i1+3);
      dy = gsl_vector_get(ee,i1+4);
      dz = gsl_vector_get(ee,i1+5);
      Node node1 = Node(x,y,z,dx,dy,dz);
      eb[i].UpdateNodes(node0,node1);
    }
  }

  void resetForce()
  {
    gsl_vector_set_zero(qe);
    gsl_vector_set_zero(qg);
    gsl_vector_set_zero(qh);
    gsl_vector_set_zero(qc);
  }
  
  void apply_bc()
  {
    // Setup the force vectors
    gsl_vector_set_zero(rhs);
    gsl_vector_sub(rhs, qe);
    gsl_vector_add(rhs, qg);
    gsl_vector_add(rhs, qh);
    gsl_vector_add(rhs, qc);
    // Set up the boundary condition
    for(int i=0; i<6; i++)
      gsl_vector_set(rhs, i, 0.0);
  }

  void readData(char *fname, char *fvname)
  {
    FILE *fp;
    fp = fopen(fname,"r+");
    if(fp==NULL){
      std::cout << "Could not open file " << fname << std::endl;
    } else{ 
      // Read data
      int nd = 6*(ne+1);
      double data;
      for(int i=0; i<nd; i++){
	fscanf(fp, "%lf\n", &data);
	gsl_vector_set(ee,i,data);
      }
    }
    fclose(fp);

    // Read velocity
    fp = fopen(fvname,"r+");
        if(fp==NULL){
      std::cout << "Could not open file " << fvname << std::endl;
    } else{ 
      // Read data
      int nd = 6*(ne+1);
      double data;
      for(int i=0; i<nd; i++){
	fscanf(fp, "%lf\n", &data);
	gsl_vector_set(u,i,data);
      }
    }
    fclose(fp);
    updateBeam();
  }

  // Advance in time for one step
  void advance(double dt, int steps)
  {
    // elasticForces();
#if 0
  FILE *fp;
  // test with given ee
  fp = fopen("ve.txt","r+");
  double e1[66];
  for(int i=0; i<66; i++){
     fscanf(fp,"%lf\n",&e1[i]);
     // std::cout << "Element ee[ " << i <<" ] = " << e1[i] << std::endl;
     gsl_vector_set(ee,i,e1[i]);
  }
  fclose(fp);
  setBeam(e1);

  elasticForces();

  fp = fopen("qev.txt","w+");
  gsl_vector_fprintf(fp,qe,"%9.6e");
  fclose(fp);

  double pf[3] = {0,0,10};
  pointLoad(9,eb[9].getLength_l(),&pf[0]);
#endif
  }

  void advanceRK2Stage1(double dt, double *xs, double **fs)
  {
    double alpha = 0.;
    double hdt;
    double pf[3];
    int ne;
    ne = getNE();
    hdt = 0.5*dt;

    // Initialize the stage
    gsl_vector_memcpy(eo,ee);
    gsl_vector_memcpy(uo,u);
    // stage 1
    resetForce();
    for(int i=1; i<ne; i++){
      pf[0] = fs[0][i];
      pf[1] = fs[1][i];
      pf[2] = fs[2][i];
      pointLoad(i,xs[i],&pf[0]);
    }
    elasticForces();
    apply_bc();
    gsl_linalg_LU_band_solve(11,11,MB,piv,rhs,sol);
    gsl_blas_daxpy(-alpha,u,sol);
    gsl_blas_daxpy(hdt,sol,u);
    gsl_blas_daxpy(hdt,u,ee);
  }

  void advanceRK2Stage2(double dt, double *xs, double **fs)
  {
    double alpha = 0.;
    double pf[3];
    int ne;
    ne  = getNE();
    resetForce();
    for(int i=1; i<ne; i++){
      pf[0] = fs[0][i];
      pf[1] = fs[1][i];
      pf[2] = fs[2][i];
      pointLoad(i,xs[i],&pf[0]);
    }
    elasticForces();
    apply_bc();
    gsl_linalg_LU_band_solve(11,11,MB,piv,rhs,sol);
    gsl_blas_daxpy(-alpha,u,sol);
    gsl_vector_memcpy(u,uo);
    gsl_blas_daxpy(dt,sol,u);
    gsl_vector_memcpy(ee,eo);
    gsl_blas_daxpy(dt,u,ee);
    // Update beam
    updateBeam();
  }

  void advanceRK2(double dt)
  {
    double alpha = 2.0;
    double hdt;
    double pf[3] = {0,0,10};
    hdt = 0.5*dt;
    // Initialize the stage
    gsl_vector_memcpy(eo,ee);
    gsl_vector_memcpy(uo,u);
    // stage 1
    resetForce();
    pointLoad(9,eb[9].getLength_l(),&pf[0]);
    elasticForces(); // use ee vector to calculate the elastic force
    apply_bc();
    gsl_linalg_LU_band_solve(11,11,MB,piv,rhs,sol); // sol = f
    gsl_blas_daxpy(-alpha,u,sol); // add the damping term
    gsl_blas_daxpy(hdt,sol,u); // u = u + 0.5*dt*f
    gsl_blas_daxpy(hdt,u,ee); //  ee = ee + 0.5*dt*u
    // stage 2
    resetForce();
    pointLoad(9,eb[9].getLength_l(),&pf[0]);
    elasticForces();
    apply_bc();
    gsl_linalg_LU_band_solve(11,11,MB,piv,rhs,sol);
    gsl_blas_daxpy(-alpha,u,sol); // add the damping term
    gsl_vector_memcpy(u,uo); 
    gsl_blas_daxpy(dt,sol,u); // update u
    gsl_vector_memcpy(ee,eo);
    gsl_blas_daxpy(dt,u,ee); // ee updated
    // update
    updateBeam();
  }

  void advanceEuler(double dt)
  {
    double alpha = 1.;
    double pf[3] = {0,0,10};
    resetForce();
    pointLoad(9,eb[9].getLength_l(), &pf[0]);
    elasticForces();
    apply_bc();
    // solver
    gsl_linalg_LU_band_solve(11,11,MB,piv,rhs,sol); // sol = f
    // update the velocity vector		
    gsl_blas_daxpy(-alpha,u,sol);
    gsl_blas_daxpy(dt,sol,u); // u = u + Dt*f
    // update the position vector
    gsl_blas_daxpy(dt,u,ee); // ee = ee + Dt*u
    // update the vector
    updateBeam();
#if 0
    FILE *fp;

    fp=fopen("u.txt","w+");
    gsl_vector_fprintf(fp,u,"%9.6e");
    fclose(fp);

    fp = fopen("sol.txt","w+");
    gsl_vector_fprintf(fp,sol,"%9.6e");
    fclose(fp);

    fp = fopen("eeu.txt","w+");
    gsl_vector_fprintf(fp,ee,"%9.6e");
    fclose(fp);
#endif
  }

  void output(char *fname, char *fvname)
  {
    FILE *fp;
    // output elements
    fp = fopen(fname,"w+");
    gsl_vector_fprintf(fp,ee,"%12.9e");
    fclose(fp);
    // output velocity
    fp = fopen(fvname,"w+");
    gsl_vector_fprintf(fp,u,"%12.9e");
    fclose(fp);
  }

  int getNE()
  {
  return ne;
  }

  gsl_vector* getEE()
  {
  return ee;
  }
 
  gsl_vector* getU()
  {
  return u;
  }

  double getLength()
  {
    return length;
  }

  double getWidth()
  {
    return width;
  }

  double getThickness()
  {
    return thickness;
  }

};

#endif
