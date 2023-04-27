// Create a beam and FSInterface

struct vec3 pos;
pos.x = 0.25;
pos.y = 0.05;
pos.z = 0.2;

Beam bm(120, 0.35, pos, 0, 0.1, 0.02, 6.67e6, 0.4, 1.e3);
FSInterface fsi(bm,6);

double **xl;
double **xc;
double **xr;
double *pl;
double *pr;
double **uc;
int ntp;

xl = fsi.getXl();
xc = fsi.getXc();
xr = fsi.getXr();
pl = fsi.getPl();
pr = fsi.getPr();
uc = fsi.getUc();
ntp = fsi.getTNP();

