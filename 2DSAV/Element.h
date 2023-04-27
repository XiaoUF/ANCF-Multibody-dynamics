#ifndef ELEMENT_H
#define ELEMENT_H

#include "node.h"

class Element{
 private:
  int index;
  Node node0;
  Node node1;
  double nu;              // Poisson's ratio
  double E;               // Young's Modulus
  double rho;             // Density of the element
  double l;               // length of the element, a thin beam element
  double w;               // width of the element
  double t;               // thickness of the element
  double m;               // mass of the element
  double collisionRadius; // This is used for collision detection
  double im;             // moment of inertia
  double area;           // area of the cross-section
  double S[4];            // Shape function
  
 public:
  // Constructor
  Element();
  Element(Node, Node);
  Element(Node, Node, double, double, double);
  Element(Node, Node, double, double);
  Element(Node, Node, double, double, double, double, double);
  // Destructor
  ~Element() {};
  
  Node getNode0(){
    return this->node0;
  }

  double getMomentInertia(){
    return this->im;
  }

  double getArea(){
    return this->area;
  }

  Node getNode1(){
    return this->node1;
  }

  double getYoungsModulus(){
    return this->E;
  }

  double getNu(){
    return this->nu;
  }

  double getWidth(){
    return this->w;
  }

  double getThickness(){
    return this->t;
  }

  double getDensity(){
    return this->rho;
  }

  double getLength_l(){
    return this->l;
  }

  double getLength(Node node1, Node node2){
    return sqrt(pow(node1.x-node2.x,2)+pow(node1.y-node2.y,2)
		+pow(node1.z-node2.z,2));
  }

  double getElementIndex(){
    return this->index;
  }

  double getCollisionRadius(){
    return this->collisionRadius;
  }

  int setRadius(double r){
    this->collisionRadius = r;
    return 0;
  }

  int setNu(double nu){
    this->nu = nu;
    return 0;
  }

  int setDensity(double density){
    this->rho = density;
    return 0;
  }

  int setElasticModulus(double E){
    this->E = E;
    return 0;
  }

  int setElementIndex(int index){
    this->index = index;
    return 0;
  }

  int setCollisionRadius(double collisionRadius){
    this->collisionRadius = collisionRadius;
    return 0;
  }

  void printInfo(){
    node0.printInfo();
    node1.printInfo();
    std::cout << "Length of the element is " << l << std::endl;
  }
  
  void UpdateNodes(Node, Node);
  double *shapeFunc(double, int);
  void shapeFunc(double, int, double *);
  void S_e(double *, double *, double *);
  void S_e(double *, double *);
  void eps_eps_e(double *, double *, double *);
  void eps_eps_e(double *, double *);
  void kappa_kappa_e(double *, double *, double, double, double *, double *, double *, double *);
};


#endif
