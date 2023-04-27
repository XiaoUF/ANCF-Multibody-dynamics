itr = itr+1;
// Update the fsi points
vector pos;
label cellID;
interpolationCellPoint<scalar> interp(p);

int ncouple = 500;
double dt = runTime.deltaT().value()/(static_cast<double> (ncouple));
double rho = 1000.;

fsi1.calc_points();
fsi1.calc_slopes();

fsi2.calc_points();
fsi2.calc_slopes();

label dummy(1);
Pstream::scatter(dummy); // act as MPI barrier

for(int i=0; i<ntp; i++){
pos = vector(xl1[0][i], xl1[1][i], xl1[2][i]);
cellID = mesh.findCell(pos);
scalar pmin;
scalar pmax;
scalar pp;
if(cellID>-1) 
  pp = rho*interp.interpolate(pos,cellID);
else
  pp = 0.0;
pmin = pp;
pmax = pp;
reduce(pmin, minOp<scalar>());
reduce(pmax, maxOp<scalar>());
pl1[i] = pmin+pmax;
}
// All reduce
for(int i=0; i<ntp; i++){
pos = vector(xr1[0][i], xr1[1][i], xr1[2][i]);
cellID = mesh.findCell(pos);
scalar pp;
scalar pmin;
scalar pmax;
if(cellID>-1)
  pp = rho*interp.interpolate(pos,cellID);
else
  pp = 0.0;
pmax = pp;
pmin = pp;
reduce(pmax,maxOp<scalar>());
reduce(pmin,minOp<scalar>());
pr1[i] = pmin+pmax;
}

// Second Beam
for(int i=0; i<ntp; i++){
pos = vector(xl2[0][i], xl2[1][i], xl2[2][i]);
cellID = mesh.findCell(pos);
scalar pmin;
scalar pmax;
scalar pp;
if(cellID>-1)
  pp = rho*interp.interpolate(pos,cellID);
else
  pp = 0.0;
pmin = pp;
pmax = pp;
reduce(pmin, minOp<scalar>());
reduce(pmax, maxOp<scalar>());
pl2[i] = pmin+pmax;
}

for(int i=0; i<ntp; i++){
pos = vector(xr2[0][i], xr2[1][i], xr2[2][i]);
cellID = mesh.findCell(pos);
scalar pp;
scalar pmin;
scalar pmax;
if(cellID>-1)
  pp = rho*interp.interpolate(pos,cellID);
else
  pp = 0.0;
pmax = pp;
pmin = pp;
reduce(pmax,maxOp<scalar>());
reduce(pmin,minOp<scalar>());
pr2[i] = pmin+pmax;
}
// Allocate space
gsl_vector *pp;
gsl_vector *pm;
gsl_movstat_workspace *w;
pp = gsl_vector_alloc(ntp);
pm = gsl_vector_alloc(ntp);
w = gsl_movstat_alloc(10);
// Left side
for(int i=0; i<ntp; i++)
   gsl_vector_set(pp,i,pl1[i]);
gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,pp,pm,w);
for(int i=0; i<ntp; i++)
   pl1[i] = gsl_vector_get(pm,i);

for(int i=0; i<ntp; i++)
   gsl_vector_set(pp,i,pl2[i]);
gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,pp,pm,w);
for(int i=0; i<ntp; i++)
   pl2[i] = gsl_vector_get(pm,i);
// Right side
for(int i=0; i<ntp; i++)
   gsl_vector_set(pp,i,pr1[i]);
gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,pp,pm,w);
for(int i=0; i<ntp; i++)
   pr1[i] = gsl_vector_get(pm,i);

for(int i=0; i<ntp; i++)
   gsl_vector_set(pp,i,pr2[i]);
gsl_movstat_mean(GSL_MOVSTAT_END_PADVALUE,pp,pm,w);
for(int i=0; i<ntp; i++)
   pr2[i] = gsl_vector_get(pm,i);
// Free gsl_vectors
gsl_vector_free(pp);
gsl_vector_free(pm);
gsl_movstat_free(w);

// Assume small movement in one timestep, and hence the hydrodynamic force
// does not change in one cycle
fsi1.hydroForce();
fsi2.hydroForce();

for(int ic=0; ic<ncouple; ic++) {
// Stage 1
 bm1.advanceRK2Stage1(dt,fsi1.getXs(),fsi1.getFs());
 bm2.advanceRK2Stage1(dt,fsi2.getXs(),fsi2.getFs());
// Stage 2
 bm1.advanceRK2Stage2(dt,fsi1.getXs(),fsi1.getFs());
 bm2.advanceRK2Stage2(dt,fsi2.getXs(),fsi2.getFs());
}

Pstream::scatter(dummy); // act as MPI barrier
//std::cout << "MPI Barrier" << std::endl;

// IO
int nout = 100;
// Test output
if(Pstream::master()){
if(itr%nout == 0){
/*
   FILE *fp;
   char fname[64];
   sprintf(fname,"Pressure_%d.txt",itr/nout);
   fp = fopen(fname,"w+");
   for(int i=0; i<ntp; i++){
     fprintf(fp,"%12.9e\t%12.9e\t%12.9e\t%12.9e\t%12.9e\n",xl[0][i],xl[1][i],xl[2][i],pl[i],pr[i]);
   }
   fclose(fp);
*/
  char fe[64];
  char fv[64];
  sprintf(fe,"Structure/eeu1_%d.txt",itr/nout);
  sprintf(fv,"Structure/uu1_%d.txt",itr/nout);
  bm1.output(fe,fv);

  sprintf(fe,"Structure/eeu2_%d.txt",itr/nout);
  sprintf(fv,"Structure/uu2_%d.txt",itr/nout);
  bm2.output(fe,fv);
}
}
// Setup the mask
double hw = 0.5*fsi1.getThickness();
forAll(mesh.C(),i)
{
  vector coords = mesh.C()[i];
  double dist1,dist2;
  double uu1[3],uu2[3];
  fsi1.MaskU(coords[0],coords[1],coords[2],dist1,&uu1[0]);
  fsi2.MaskU(coords[0],coords[1],coords[2],dist2,&uu2[0]);
  scalar vsx, vsz;
  D[i] = 0.5*K.value()*(1.0-Foam::tanh(2000.*(dist1-hw)))
       + 0.5*K.value()*(1.0-Foam::tanh(2000.*(dist2-hw)));
  vsx  = 0.5*uu1[0]*(1.0-Foam::tanh(2000.*(dist1-hw)))
       + 0.5*uu2[0]*(1.0-Foam::tanh(2000.*(dist2-hw)));
  vsz  = 0.5*uu1[2]*(1.0-Foam::tanh(2000.*(dist1-hw)))
       + 0.5*uu2[2]*(1.0-Foam::tanh(2000.*(dist2-hw)));
  Uf[i] = vector(vsx, 0.0, vsz);
 // Uf[i] = vector(0.,0.,0.);
}

Pstream::scatter(dummy); // act as MPI barrier

