
// Read Darcy's coefficient
Info << "\nReading Darcy's coefficient"<< endl;
IOdictionary DarcyProperties
(
 IOobject
 (
  "DarcyProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );

const dimensionedScalar K(DarcyProperties.lookup("K"));
//const dimensionedScalar XC(DarcyProperties.lookup("xc"));
//const dimensionedScalar YC(DarcyProperties.lookup("yc"));
//const dimensionedScalar Radius(DarcyProperties.lookup("rc"));
//const dimensionedVector f_cur(DarcyProperties.lookup("fc"));
// This is used to create the coefficients required for Porous media
volScalarField D
(
 IOobject
 (
  "D_epsilon",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("0",K.dimensions(),0.0)
 );

volVectorField Uf
(
 IOobject
 (
  "Uf",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedVector("0",U.dimensions(),vector(0.0,0.0,0.0))
 );


//create the permeability 
Info << "Initialize the permeability...\n" << endl;
forAll(mesh.C(),i)
{
//  scalar dist = Foam::sqrt(a*a+b*b);
  Uf[i] = vector(0.0, 0.0, 0.0);
  D[i]  = K.value();
}

