itr = itr+1;

double DX = 0.0005;
double DZ = 0.0005;
double Pi = Foam::constant::mathematical::pi;

label dummy(1);
Pstream::scatter(dummy); // act as MPI barrier

// Update the fsi points
// Initialize the data
for(int i=0; i<TNB; i++)
  for(int j=0; j<NBP; j++){
    xbc[i][j] = 0.0;
    ybc[i][j] = 0.0;
    zbc[i][j] = 0.0;
    ubc[i][j] = 0.0;
    vbc[i][j] = 0.0;
    wbc[i][j] = 0.0;
    xlc[i][j] = 0.0;
    ylc[i][j] = 0.0;
    zlc[i][j] = 0.0;
    xrc[i][j] = 0.0;
    yrc[i][j] = 0.0;
    zrc[i][j] = 0.0;    
    plc[i][j] = 0.0;
    prc[i][j] = 0.0;
  }

Pstream::scatter(dummy); // act as MPI barrier
// Loop over each beam
int ntp = fsi[0].getTNP();;
int ncouple = 100;
double dt = runTime.deltaT().value()/(static_cast<double> (ncouple));
double rho = 1000.;

for(int ib=0; ib<TNB; ib++)
  for(int ip=0; ip<NP; ip++){
    xp[ib][ip] = 0.0;
    zp[ib][ip] = 0.0;
    up[ib][ip] = 0.0;
    wp[ib][ip] = 0.0;
}
// Setup particle locations and velocities
for(int ib=0; ib<nb; ib++){
  int bid;
  bid = nb*pid+ib;
  fsi[ib].GetParticles(&xp[bid][0],&zp[bid][0],&up[bid][0],&wp[bid][0],NPE);
}
// Setup AABB box
for(int ib=0; ib<nb; ib++){
  int bid;
  bid = nb*pid+ib;
  for(int ie=0; ie<ne; ie++){
    double xlo,xhi,zlo,zhi;
    fsi[ib].GetBox(ie,xlo,zlo,xhi,zhi);
    int idx;
    idx = bid*ne+ie;
    box[idx][0] = xlo;
    box[idx][1] = zlo;
    box[idx][2] = xhi;
    box[idx][3] = zhi; 
  }
}

for(int ib=0; ib<nb; ib++){
  double **xl;
  double **xr;

  int bid;
  
  bid = nb*pid+ib;
  xl = fsi[ib].getXl();
  xr = fsi[ib].getXr();

  fsi[ib].calc_points();
  fsi[ib].calc_slopes();
  // Set up data for the communication
  for(int i=0; i<ntp; i++){
    xlc[bid][i] = xl[0][i];
    ylc[bid][i] = xl[1][i];
    zlc[bid][i] = xl[2][i];
    xrc[bid][i] = xr[0][i];
    yrc[bid][i] = xr[1][i];
    zrc[bid][i] = xr[2][i];
  }
 }

int count = TNB*NBP;
MPI_Allreduce(MPI_IN_PLACE,xlc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,ylc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,zlc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,xrc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,yrc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,zrc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

count = TNB*NP;
MPI_Allreduce(MPI_IN_PLACE,xp[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,zp[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,up[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,wp[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

count = TNB*ne*4;
MPI_Allreduce(MPI_IN_PLACE,box[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

//if(Pstream::master())

// Communication to get the pressure information
for(int ib=0; ib<TNB; ib++) {
  for(int i=0; i<ntp; i++){
    // Get  pressure for left side
    scalar pp;
    pp = 0;
    forAll(mesh.C(),j)
      {
	vector coords = mesh.C()[j];
	double xx,zz;
	double deltax, deltaz;
	xx = std::fabs(coords[0]-xlc[ib][i])/DX;
	zz = std::fabs(coords[2]-zlc[ib][i])/DZ;
	deltax = 0.0;
	deltaz = 0.0;
	if(xx<=2.5 && zz<=2.5){
	  // Get deltax
	  if(xx<=0.5)
	    deltax = 0.375+Pi/32-0.25*xx*xx;
	  else if(xx<=1.5)
	    deltax = 0.25+0.125*(1.-xx)*std::sqrt(-2.+8.*xx-4.*xx*xx)-0.125*std::asin(std::sqrt(2)*(xx-1.));
	  else if(xx<=2.5)
	    deltax = 17./16.-Pi/64.-0.75*xx+0.125*xx*xx+(xx-2.)*std::sqrt(16*xx-4.*xx*xx-14.)/16.+std::asin(std::sqrt(2)*(xx-2.))/16.;
	  else
	    deltax = 0.0;
	  // Get deltaz
	  if(zz<=0.5)
	    deltaz = 0.375+Pi/32-0.25*zz*zz;
	  else if(zz<=1.5)
	    deltaz = 0.25+0.125*(1.-zz)*std::sqrt(-2.+8.*zz-4.*zz*zz)-0.125*std::asin(std::sqrt(2)*(zz-1.));
	  else if(zz<=2.5)
	    deltaz = 17./16.-Pi/64.-0.75*zz+0.125*zz*zz+(zz-2.)*std::sqrt(16*zz-4.*zz*zz-14.)/16.+std::asin(std::sqrt(2)*(zz-2.))/16.;
	  else
	    deltaz = 0.0;
	}
	pp = pp+deltax*deltaz*p[j]; 
      }
    reduce(pp,sumOp<scalar>());
    //std::cout << "pp = " << pp << std::endl;
    plc[ib][i]=rho*pp;
  }

  for(int i=0; i<ntp; i++){
    // use the regularized delta function
    scalar pp;
    pp = 0;
    forAll(mesh.C(),j)
      {
	vector coords = mesh.C()[j];
	double xx,zz;
	double deltax, deltaz;
	xx = std::fabs(coords[0]-xrc[ib][i])/DX;
	zz = std::fabs(coords[2]-zrc[ib][i])/DZ;
	deltax = 0.0;
	deltaz = 0.0;
	if(xx<=2.5 && zz<=2.5){
	  // Get deltax
	  if(xx<=0.5)
	    deltax = 0.375+Pi/32-0.25*xx*xx;
	  else if(xx<=1.5)
	    deltax = 0.25+0.125*(1.-xx)*std::sqrt(-2.+8.*xx-4.*xx*xx)-0.125*std::asin(std::sqrt(2)*(xx-1.));
	  else if(xx<=2.5)
	    deltax = 17./16.-Pi/64.-0.75*xx+0.125*xx*xx+(xx-2.)*std::sqrt(16*xx-4.*xx*xx-14.)/16.+std::asin(std::sqrt(2)*(xx-2.))/16.;
	  else
	    deltax = 0.0;
	  // Get deltaz
	  if(zz<=0.5)
	    deltaz = 0.375+Pi/32-0.25*zz*zz;
	  else if(zz<=1.5)
	    deltaz = 0.25+0.125*(1.-zz)*std::sqrt(-2.+8.*zz-4.*zz*zz)-0.125*std::asin(std::sqrt(2)*(zz-1.));
	  else if(zz<=2.5)
	    deltaz = 17./16.-Pi/64.-0.75*zz+0.125*zz*zz+(zz-2.)*std::sqrt(16*zz-4.*zz*zz-14.)/16.+std::asin(std::sqrt(2)*(zz-2.))/16.;
	  else
	    deltaz = 0.0;
	}
	pp = pp+deltax*deltaz*p[j]; 
      }
    reduce(pp,sumOp<scalar>());
    //std::cout << "pr = " << pp << std::endl;
    prc[ib][i]=rho*pp;
  }
 }

Pstream::scatter(dummy); // act as MPI barrier

if(Pstream::master())
  std::cout << "Advance the beam using hydro-force " << std::endl;
// Assume small movement in one timestep, and hence the hydrodynamic force
// does not change in one cycle
for(int ib=0; ib<nb; ib++){
  double *pr;
  double *pl;
  double **xc;
  double **uc;
  int bid;
  bid = nb*pid+ib;
  pl = fsi[ib].getPl();
  pr = fsi[ib].getPr();
  xc = fsi[ib].getXc();
  uc = fsi[ib].getUc();

  for(int i=0; i<ntp; i++){
    pl[i] = plc[bid][i];
    pr[i] = prc[bid][i];
  }
  fsi[ib].hydroForce();

  for(int ic=0; ic<ncouple; ic++) {
    // Stage 1
    ba[ib].advanceRK2Stage1(dt,fsi[ib].getXs(),fsi[ib].getFs());
    // Stage 2
    ba[ib].advanceRK2Stage2(dt,fsi[ib].getXs(),fsi[ib].getFs());
  }
  // update beam center
  fsi[ib].calc_points();
  for(int i=0; i<ntp; i++){
    xbc[bid][i] = xc[0][i];
    ybc[bid][i] = xc[1][i];
    zbc[bid][i] = xc[2][i];
    ubc[bid][i] = uc[0][i];
    vbc[bid][i] = uc[1][i];
    wbc[bid][i] = uc[2][i];
  }
  // IO
  int nout = 1000;
  if(itr%nout == 0){
    char fe[64];
    char fv[64];
    sprintf(fe,"Structure/eeu%d_%d.txt",bid,itr/nout);
    sprintf(fv,"Structure/uu%d_%d.txt",bid,itr/nout);
    ba[ib].output(fe,fv);
  }
 }

Pstream::scatter(dummy); // act as MPI barrier
MPI_Allreduce(MPI_IN_PLACE,xbc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,ybc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,zbc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,ubc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,vbc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE,wbc[0],count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

Pstream::scatter(dummy); // act as MPI barrier

// Setup the mask
double hw = 0.5*thickness;
forAll(mesh.C(),i)
{
  vector coords = mesh.C()[i];
  double dist;
  double uu[3];
  //    Get the distance
  dist = 1.e6;
  scalar vsx,vsz;
  for(int j=0; j<TNB; j++) {
    for(int k=0; k<ntp; k++){
      double dx,dy,dz;
      double d2;
      dx = coords[0]-xbc[j][k];
      dy = coords[1]-ybc[j][k];
      dz = coords[2]-zbc[j][k];
      d2 = std::sqrt(dx*dx+dz*dz);
      if(d2<dist){
	dist = d2;
	uu[0] = ubc[j][k];
	uu[1] = vbc[j][k];
	uu[2] = wbc[j][k];
      }
    }
  }
  D[i] = 0.5*K.value()*(1.0-Foam::tanh(3000.*(dist-hw)));
  vsx  = 0.5*uu[0]*(1.0-Foam::tanh(3000.*(dist-hw)));
  vsz  = 0.5*uu[2]*(1.0-Foam::tanh(3000.*(dist-hw)));
  Uf[i] = vector(vsx, 0.0, vsz);

}

Pstream::scatter(dummy); // act as MPI barrier
