#include "Element.h"

Element::Element(){
    // Default constructor
    this->node0 = Node(0,0,0,1,0,0);
    this->node1 = Node(1,0,0,1,0,0);
    // geometric measures of the vegetation blade element
    this->l   = 1.0;
    this->w   = 0.1;
    this->t   = 0.02;
    // Mechanical properties of the vegetation blade
    this->E   = 2.e7;
    this->nu  = 0.3;
    this->rho = 1150.0;
    this->m   = rho*l*w*t;
    this->im  = w*t*t*t/12.0;
    this->area = w*t;
    collisionRadius = 0.0; 
  }
Element::Element(Node node0, Node node1)
    {
      // Constructor from two nodes
      this->node0 = node0;
      this->node1 = node1;
      this->l = getLength(node0,node1);
      this->nu  = 0.3;
      this->E   = 2.e7;
      this->rho = 1150.0;
      this->w   = 0.1;
      this->t   = 0.02;
      this->m   = rho*l*w*t;
      this->im  = w*t*t*t/12.0;
      this->area = w*t;
      collisionRadius = 0.0;    
    }
Element::Element(Node node0, Node node1, double E, double nu, double rho)
    {
      // constructor from two nodes, and material properties
      this->node0 = node0;
      this->node1 = node1;
      this->l     = getLength(node0,node1);
      this->nu    = nu;
      this->E     = E;
      this->rho   = rho;
      this->w   = 0.1;
      this->t   = 0.02;
      this->m   = rho*l*w*t;
      this->im  = w*t*t*t/12.0;
      this->area = w*t;
      collisionRadius = 0.0;   
    }
  
Element::Element(Node node0, Node node1, double w, double t)
    {
      // constructor from two nodes and vegetation blade measures
      this->node0 = node0;
      this->node1 = node1;
      this->l     = getLength(node0,node1);
      this->E     = 2.e7;
      this->nu    = 0.3;
      this->w     = w;
      this->t     = t;
      this->rho   = 1150.0;
      this->m     = rho*l*w*t;
      this->im    = w*t*t*t/12.0;
      this->area = w*t;
      collisionRadius = 0.0;
    }
  
Element::Element(Node node0, Node node1, double w, double t, double E, double nu, double rho)
    {
      // constructor 
      this->node0 = node0;
      this->node1 = node1;
      this->l     = getLength(node0,node1);
      this->nu    = nu;
      this->E     = E;
      this->rho   = rho;
      this->w     = w;
      this->t     = t;
      this->m     = rho*l*w*t;
      this->im    = w*t*t*t/12.0;
      this->area = w*t;
      collisionRadius = 0.0;
    }

void Element::UpdateNodes(Node node0, Node node1)
{
  this->node0.x  = node0.x;
  this->node0.y  = node0.y;
  this->node0.z  = node0.z;
  this->node0.dx = node0.dx;
  this->node0.dy = node0.dy;
  this->node0.dz = node0.dz;

  this->node1.x  = node1.x;
  this->node1.y  = node1.y;
  this->node1.z  = node1.z;
  this->node1.dx = node1.dx;
  this->node1.dy = node1.dy;
  this->node1.dz = node1.dz;

}

  
double * Element::shapeFunc(double x,int order)
  {
    double xi = x/l;
    double li  = 1./l;
    double li2 = li*li;
    if(order < 0 || order > 2){
      std::cout << "Order should be an integer in the range of [0, 2]"
		<< std::endl;
      return NULL;
    }
    // std::cout << "xi = " << xi <<std::endl;
    switch(order){
    case 0:
      S[0] = 1.0-3.*xi*xi + 2.*xi*xi*xi;
      S[1] = l*(xi - 2.*xi*xi + xi*xi*xi);
      S[2] = 3.*xi*xi - 2.*xi*xi*xi;
      S[3] = l*(-xi*xi + xi*xi*xi);
      break;
    case 1:
      S[0] = (6.*xi*xi - 6*xi)*li;
      S[1] = 1. - 4.*xi + 3.*xi*xi;
      S[2] = (-6.*xi*xi + 6*xi)*li;
      S[3] = -2.*xi + 3.*xi*xi;
      break;
    case 2:
      S[0] = (12.*xi - 6.)*li2;
      S[1] = (-4. + 6.*xi)*li;
      S[2] = (6. - 12.*xi)*li2;
      S[3] = (-2. + 6.*xi)*li;
      break;
    default :
      S[0] = 0.;
      S[1] = 0.;
      S[2] = 0.;
      S[3] = 0.;
    }
    return &S[0];
  }

void Element::shapeFunc(double x,int order, double *Sx)
  {
    double xi = x/l;
    double li  = 1./l;
    double li2 = li*li;
    if(order < 0 || order > 2){
      std::cout << "Order should be an integer in the range of [0, 2]"
		<< std::endl;
    }
    // std::cout << "xi = " << xi <<std::endl;
    switch(order){
    case 0:
      Sx[0] = 1.0-3.*xi*xi + 2.*xi*xi*xi;
      Sx[1] = l*(xi - 2.*xi*xi + xi*xi*xi);
      Sx[2] = 3.*xi*xi - 2.*xi*xi*xi;
      Sx[3] = l*(-xi*xi + xi*xi*xi);
      break;
    case 1:
      Sx[0] = (6.*xi*xi - 6*xi)*li;
      Sx[1] = 1. - 4.*xi + 3.*xi*xi;
      Sx[2] = (-6.*xi*xi + 6*xi)*li;
      Sx[3] = -2.*xi + 3.*xi*xi;
      break;
    case 2:
      Sx[0] = (12.*xi - 6.)*li2;
      Sx[1] = (-4. + 6.*xi)*li;
      Sx[2] = (6. - 12.*xi)*li2;
      Sx[3] = (-2. + 6.*xi)*li;
      break;
    default :
      Sx[0] = 0.;
      Sx[1] = 0.;
      Sx[2] = 0.;
      Sx[3] = 0.;
    }
  }


void Element::eps_eps_e(double *Sx, double *res)
{
  double e[12] = {node0.x,node0.y,node0.z,node0.dx,node0.dy,node0.dz,node1.x,node1.y,node1.z,node1.dx,node1.dy,node1.dz};
  double t2 = Sx[0];
  double t3 = e[0];
  double t4 = t2*t3;
  double t5 = Sx[1];
  double t6 = e[3];
  double t7 = t5*t6;
  double t8 = Sx[2];
  double t9 = e[6];
  double t10 = t8*t9;
  double t11 = Sx[3];
  double t12 = e[9];
  double t13 = t11*t12;
  double t14 = t4 + t7 + t10 + t13;
  double t17 = e[1];
  double t18 = t2*t17;
  double t19 = e[4];
  double t20 = t5*t19;
  double t21 = e[7];
  double t22 = t8*t21;
  double t23 = e[10];
  double t24 = t11*t23;
  double t15 = t18 + t20 + t22 + t24;
  double t29 = e[2];
  double t30 = t2*t29;
  double t31 = e[5];
  double t32 = t5*t31;
  double t33 = e[8];
  double t34 = t8*t33;
  double t35 = e[11];
  double t36 = t11*t35;
  double t16 = t30 + t32 + t34 + t36;
  double t25 = t14*t14;
  double t26 = t25*0.5;
  double t27 = t15*t15;
  double t28 = 0.5*t27;
  double t37 = t16*t16;
  double t38 = t37*0.5;
  double t39 = t26 + t28 + t38 -0.5;

  res[0] = t2*t14*t39;
  res[1] = t2*t15*t39;
  res[2] = t2*t16*t39;
  res[3] = t5*t14*t39;
  res[4] = t5*t15*t39;
  res[5] = t5*t16*t39;
  res[6] = t8*t14*t39;
  res[7] = t8*t15*t39;
  res[8] = t8*t16*t39;
  res[9] = t11*t14*t39;
  res[10] = t11*t15*t39;
  res[11] = t11*t16*t39;
}

void Element::eps_eps_e(double *Sx, double *e, double *res)
{
  double t2 = Sx[0];
  double t3 = e[0];
  double t4 = t2*t3;
  double t5 = Sx[1];
  double t6 = e[3];
  double t7 = t5*t6;
  double t8 = Sx[2];
  double t9 = e[6];
  double t10 = t8*t9;
  double t11 = Sx[3];
  double t12 = e[9];
  double t13 = t11*t12;
  double t14 = t4 + t7 + t10 + t13;
  double t17 = e[1];
  double t18 = t2*t17;
  double t19 = e[4];
  double t20 = t5*t19;
  double t21 = e[7];
  double t22 = t8*t21;
  double t23 = e[10];
  double t24 = t11*t23;
  double t15 = t18 + t20 + t22 + t24;
  double t29 = e[2];
  double t30 = t2*t29;
  double t31 = e[5];
  double t32 = t5*t31;
  double t33 = e[8];
  double t34 = t8*t33;
  double t35 = e[11];
  double t36 = t11*t35;
  double t16 = t30 + t32 + t34 + t36;
  double t25 = t14*t14;
  double t26 = t25*0.5;
  double t27 = t15*t15;
  double t28 = 0.5*t27;
  double t37 = t16*t16;
  double t38 = t37*0.5;
  double t39 = t26 + t28 + t38 -0.5;

  res[0] = t2*t14*t39;
  res[1] = t2*t15*t39;
  res[2] = t2*t16*t39;
  res[3] = t5*t14*t39;
  res[4] = t5*t15*t39;
  res[5] = t5*t16*t39;
  res[6] = t8*t14*t39;
  res[7] = t8*t15*t39;
  res[8] = t8*t16*t39;
  res[9] = t11*t14*t39;
  res[10] = t11*t15*t39;
  res[11] = t11*t16*t39;
}
 
void Element::kappa_kappa_e(double *Sx, double *Sxx, double coef, double inv_rx_rx_cubed, double *rx, double *rxx, double *v, double *res)
{
  double t2 = Sx[0];
  double t3 = Sxx[0];
  double t4 = rx[0];
  double t5 = v[2];
  double t6 = rxx[2];
  double t7 = t2*t6;
  double t8 = rx[2];
  double t9 = t7 - t3*t8;
  double t10 = rx[1];
  double t11 = rxx[0];
  double t12 = t3*t4;
  double t13 = v[1];
  double t14 = rxx[1];
  double t15 = t3*t10;
  double t16 = v[0];
  double t17 = t15 - t2*t14;
  double t18 = t12 - t2*t11;
  double t19 = Sx[1];
  double t20 = Sxx[1];
  double t21 = t6*t19;
  double t22 = t21 - t8*t20;
  double t23 = t4*t20;
  double t24 = t10*t20;
  double t25 = t24 - t14*t19;
  double t26 = t23 - t11*t19;
  double t27 = Sx[2];
  double t28 = Sxx[2];
  double t29 = t6*t27;
  double t30 = t29 - t8*t28;
  double t31 = t4*t28;
  double t32 = t10*t28;
  double t33 = t32 - t14*t27;
  double t34 = t31 - t11*t27;
  double t35 = Sx[3];
  double t36 = Sxx[3];
  double t37 = t6*t35;
  double t38 = t37 - t8*t36;
  double t39 = t4*t36;
  double t40 = t10*t36;
  double t41 = t40 - t14*t35;
  double t42 = t39 - t11*t35;

  res[0] = -inv_rx_rx_cubed*(t5*t17+t9*t13+coef*t2*t4);
  res[1] = inv_rx_rx_cubed*(t5*t18+t9*t16-coef*t2*t10);
  res[2] = -inv_rx_rx_cubed*(t13*t18-t16*t17+coef*t2*t8);
  res[3] = -inv_rx_rx_cubed*(t5*t25+t13*t22+coef*t4*t19);
  res[4] = inv_rx_rx_cubed*(t5*t26+t16*t22-coef*t10*t19);
  res[5] = -inv_rx_rx_cubed*(t13*t26-t16*t25+coef*t8*t19);
  res[6] = -inv_rx_rx_cubed*(t5*t33+t13*t30+coef*t4*t27);
  res[7] = inv_rx_rx_cubed*(t5*t34+t16*t30-coef*t10*t27);
  res[8] = -inv_rx_rx_cubed*(t13*t34-t16*t33+coef*t8*t27);
  res[9] = -inv_rx_rx_cubed*(t5*t41+t13*t38+coef*t4*t35);
  res[10] = inv_rx_rx_cubed*(t5*t42+t16*t38-coef*t10*t35);
  res[11] = -inv_rx_rx_cubed*(t13*t42-t16*t41+coef*t8*t35);
}

void Element::S_e(double *S, double *e, double *res)
{
  for(int i=0; i<3; i++){
    res[i] = S[0]*e[i] + S[1]*e[i+3] + S[2]*e[i+6] + S[3]*e[i+9];
  }
}

void Element::S_e(double *S, double *res)
{
  double e[12] = {node0.x,node0.y,node0.z,node0.dx,node0.dy,node0.dz,node1.x,node1.y,node1.z,node1.dx,node1.dy,node1.dz};
  for(int i=0; i<3; i++){
    res[i] = S[0]*e[i] + S[1]*e[i+3] + S[2]*e[i+6] + S[3]*e[i+9];
  }
}

