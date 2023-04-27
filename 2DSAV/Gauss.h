#ifndef GAUSS_H
#define GAUSS_H

#include <math.h>

class GaussPoint
{
 public:
  double x2[2] = {-sqrt(1./3.), sqrt(1./3.)};
  double w2[2] = {1, 1};
  double x3[3] = {-sqrt(3./5.), 0, sqrt(3./5.)};
  double w3[3] = {5./9., 8./9., 5./9.};
  double x4[4] = {-sqrt(3./7.+sqrt(120.)/35.),-sqrt(3./7.-sqrt(120.)/35.),
		  sqrt(3./7.-sqrt(120.)/35.), sqrt(3./7.+sqrt(120.)/35.)};
  double w4[4] = {0.5-5./(3.*sqrt(120.)), 0.5+5./(3.*sqrt(120.)),
		  0.5+5./(3.*sqrt(120.)), 0.5-5./(3.*sqrt(120.))}; 
  double x5[5] = {-sqrt(5.+2.*sqrt(10./7.))/3.,-sqrt(5.-2.*sqrt(10./7.))/3.,
		  0.,(sqrt(5.-2.*sqrt(10./7.)))/3.,sqrt(5.+2.*sqrt(10./7.))/3.};
  double w5[5] = {(322.-13.*sqrt(70.))/900.,(322.+13.*sqrt(70.))/900.,128./225.,
		  (322.+13.*sqrt(70.))/900.,(322.-13.*sqrt(70.))/900.};
};
#endif
