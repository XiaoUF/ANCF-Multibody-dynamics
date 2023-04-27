// Create a beam and FSInterface

struct vec3 pos;
pos.x = 0.15;
pos.y = 0.005;
pos.z = 0.0;

Beam bm1(80, 0.1, pos, 2, 0.01, 0.004, 2.e5, 0.3, 2.e3);
FSInterface fsi1(bm1,2);

pos.x = 0.165;
pos.y = 0.005;
pos.z = 0.0;

Beam bm2(80, 0.1, pos, 2, 0.01, 0.004, 2.e5, 0.3, 2.e3);
FSInterface fsi2(bm2,2);

double **xl1, **xl2;
double **xr1, **xr2;
double *pl1, *pl2;
double *pr1, *pr2;
int ntp;

xl1 = fsi1.getXl();
xr1 = fsi1.getXr();
pl1 = fsi1.getPl();
pr1 = fsi1.getPr();
ntp = fsi1.getTNP();

xl2 = fsi2.getXl();
xr2 = fsi2.getXr();
pl2 = fsi2.getPl();
pr2 = fsi2.getPr();
