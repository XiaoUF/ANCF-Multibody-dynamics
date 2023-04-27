#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <stdlib.h>

#define PI 3.141592653589793238

class Particle{
 private:
  int    np;   // number of particles
  double r;    // particle radius
  double E;    // Young's Modulus
  double nu;   // Poisson's ratio 
  double e;    // coefficient of restitution
  double rho;  // density
  double m;    // mass of the particle
  double Kn;   // Kn is the spring coefficient for normal force
  double gn;   // gamma_n is the damping coefficient for normal force
  double Kt;   // Kt is the spring coefficient for tangential force
  double gt;   // gamma_t is the damping coefficient for tangential force
  double **pos;  // position vector 
  double **vel;  // velocity vector
  double **fc;   // contact force vector
 public:
  // Constructor
  Particle(int np, double r){
    this->np = np;
    this->r  = r;
    // Set the particle properties
    this->E   = 2.0e4;
    this->nu  = 0.3;
    this->e   = 0.5;
    this->rho = 1150.0;
    this->m   = rho*PI*pow(r,3)*2.0/3.0;
    // Calculate the DEM coefficient
    this->Kn = 100.;
    this->gn = 10.0;
    this->Kt = 30.0;
    this->gt = 10.0;
    // allocate array
    pos = create_2d_double_array(np,3,"Position");
    vel = create_2d_double_array(np,3,"Velocity");
    fc  = create_2d_double_array(np,3,"Force");
  }
  

  Particle(int np, double r, double Kn, double gn, double Kt, double gt){
    this->np = np;
    this->r  = r;
    // Set the particle properties
    this->E   = 2.0e4;
    this->nu  = 0.3;
    this->e   = 0.5;
    this->rho = 1150.0;
    this->m   = rho*PI*pow(r,3)*2.0/3.0;
    // Calculate the DEM coefficient
    this->Kn = Kn;
    this->gn = gn;
    this->Kt = Kt;
    this->gt = gt;
    // allocate array
    pos = create_2d_double_array(np,3,"Position");
    vel = create_2d_double_array(np,3,"Velocity");
    fc  = create_2d_double_array(np,3,"Force");
  }
  // Destructor
  ~Particle() {
    // free the resources
    destroy_2d_double_array(pos);
    destroy_2d_double_array(vel);
    destroy_2d_double_array(fc);
  }

  // Memory management
  double **create_2d_double_array(int n1, int n2, const char *name)
  {
    double *data = (double *) malloc(n1*n2*sizeof(double));
    double **array = (double **) malloc(n1*sizeof(double *));

    int n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }

  void destroy_2d_double_array(double **array)
  {
    if (array == NULL) return;
    free(array[0]);
    free(array);
  }
  
  // Set propoerties of Particle
  int setKn(double Kn){
    this->Kn = Kn;
    return 0;
  }

  int setGn(double gn){
    this->gn = gn;
    return 0;
  }

  int setKt(double Kt){
    this->Kt = Kt;
    return 0;
  }

  int setGt(double gt){
    this->gt = gt;
    return 0;
  }
  
  // Interface to the structure solvers
  /* void set_particles(int nb, Beam* bm){ */
  /*   // Need both position and velocities of the beam */
  /*   // Use linear element to setup particles */
  /* } */

  void calc_CollisionForce(){

  }
  
  double **getCollisionForce(){
    return this->fc;
  }

  
};

#endif
