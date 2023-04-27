// Create a beam and FSInterface

struct vec3 pos;
pos.x = 0.15;
pos.y = 0.005;
pos.z = 0.0;

int TNB = 32; // Total number of beams
label pid = Pstream::myProcNo(); // my processor id
label n = Pstream::nProcs();     // total number of processors
int nb = TNB/n;                  // number of local beams
TNB = nb*n;

// Create workspace for communication
int NBP = 200; // make sure to hold enough space
// Workspace for mask
double **xbc;
double **ybc;
double **zbc;
double **ubc;
double **vbc;
double **wbc;
// Workspace for hydrodynamic pressure force
double **xlc;
double **ylc;
double **zlc;
double **xrc;
double **yrc;
double **zrc;
double **plc;
double **prc;

// total number of elements in 2D array
int nd = NBP*TNB;

xbc = (double **) malloc(TNB*sizeof(double *));
ybc = (double **) malloc(TNB*sizeof(double *));
zbc = (double **) malloc(TNB*sizeof(double *));
ubc = (double **) malloc(TNB*sizeof(double *));
vbc = (double **) malloc(TNB*sizeof(double *));
wbc = (double **) malloc(TNB*sizeof(double *));

xlc = (double **) malloc(TNB*sizeof(double *));
ylc = (double **) malloc(TNB*sizeof(double *));
zlc = (double **) malloc(TNB*sizeof(double *));
xrc = (double **) malloc(TNB*sizeof(double *));
yrc = (double **) malloc(TNB*sizeof(double *));
zrc = (double **) malloc(TNB*sizeof(double *));

plc = (double **) malloc(TNB*sizeof(double *));
prc = (double **) malloc(TNB*sizeof(double *));


double *xdata,*ydata,*zdata,*udata,*vdata,*wdata;
xdata = (double *) malloc(sizeof(double)*nd);
ydata = (double *) malloc(sizeof(double)*nd);
zdata = (double *) malloc(sizeof(double)*nd);
udata = (double *) malloc(sizeof(double)*nd);
vdata = (double *) malloc(sizeof(double)*nd);
wdata = (double *) malloc(sizeof(double)*nd);
int n1=0;
for (int i=0; i<TNB; i++){
  xbc[i] = &xdata[n1];
  ybc[i] = &ydata[n1];
  zbc[i] = &zdata[n1];
  ubc[i] = &udata[n1];
  vbc[i] = &vdata[n1];
  wbc[i] = &wdata[n1];
  n1 = n1+NBP;
}

double *xldata,*yldata,*zldata,*xrdata,*yrdata,*zrdata;
xldata = (double *) malloc(sizeof(double)*nd);
yldata = (double *) malloc(sizeof(double)*nd);
zldata = (double *) malloc(sizeof(double)*nd);
xrdata = (double *) malloc(sizeof(double)*nd);
yrdata = (double *) malloc(sizeof(double)*nd);
zrdata = (double *) malloc(sizeof(double)*nd);
n1=0;
for (int i=0; i<TNB; i++){
  xlc[i] = &xldata[n1];
  ylc[i] = &yldata[n1];
  zlc[i] = &zldata[n1];
  xrc[i] = &xrdata[n1];
  yrc[i] = &yrdata[n1];
  zrc[i] = &zrdata[n1];
  n1 = n1+NBP;
}

double *pldata,*prdata;
pldata = (double *) malloc(sizeof(double)*nd);
prdata = (double *) malloc(sizeof(double)*nd);
n1=0;
for (int i=0; i<TNB; i++){
  plc[i] = &pldata[n1];
  prc[i] = &prdata[n1];
  n1 = n1+NBP;
}

// Properties of the beam
double thickness = 0.002;
double width     = 0.01;
double length    = 0.1;
int  ne = 100;
double modulus = 8.e5;
double rho     = 1.e3;

//Beam bm(ne,length,pos,2,width,thickness,modulus,0.3,rho);
Beam *ba = new Beam[nb];
//int istep = 200;
for(int i=0; i<nb; i++){
  int ib;
  ib = nb*pid+i;
  pos.x = 0.15+0.015*ib;
  pos.y = 0.005;
  pos.z = 0.0;
  ba[i].SetupBeam(ne,length,pos,2,width,thickness,modulus,0.3,rho);
// Read data
/*
  char fname[32];
  char fvname[32];
  sprintf(fname,"Structure/eeu%d_%d.txt",ib,istep);
  sprintf(fvname,"Structure/uu%d_%d.txt",ib,istep);
  ba[i].readData(fname,fvname);
*/
}

// Setup the interface
FSInterface *fsi = new FSInterface[nb];
for(int i=0; i<nb; i++){
  fsi[i].setup(&ba[i],2);
}

// Setup the particles
int NPE,NP;
NPE = ceil(length/thickness/ne); // Make sure to hold enough space
NP = NPE*ne;
std::cout << "Total number of particles is " << NP <<std::endl;
double **xp;
double **zp;
double **up;
double **wp;
double **box;

int nnp = TNB*NP;
xp = (double **) malloc(TNB*sizeof(double *));
zp = (double **) malloc(TNB*sizeof(double *));
up = (double **) malloc(TNB*sizeof(double *));
wp = (double **) malloc(TNB*sizeof(double *));

double *xpdata, *zpdata, *updata, *wpdata;
xpdata = (double *) malloc(sizeof(double)*nnp);
zpdata = (double *) malloc(sizeof(double)*nnp);
updata = (double *) malloc(sizeof(double)*nnp);
wpdata = (double *) malloc(sizeof(double)*nnp);

int n3=0;
for(int i=0; i<TNB; i++){
  xp[i] = &xpdata[n3];
  zp[i] = &zpdata[n3];
  up[i] = &updata[n3];
  wp[i] = &wpdata[n3];
  n3 = n3+NP;
 }

// Setup AABB Box

int NB = ne*4;
int nnb = TNB*NB;
box = (double **) malloc(TNB*sizeof(double *));
double *boxdata;
boxdata = (double *) malloc(sizeof(double)*nnb);
n3 = 0;
for(int i=0; i<TNB; i++){
  box[i] = &boxdata[3];
  n3 = n3+NB;
 }
