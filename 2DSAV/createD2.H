
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
const dimensionedScalar KC(DarcyProperties.lookup("KC"));
const dimensionedScalar usrc(DarcyProperties.lookup("U_src"));
const dimensionedScalar xc(DarcyProperties.lookup("xc"));
const dimensionedScalar yc(DarcyProperties.lookup("yc"));

// This is used to create the coefficients required for Porous media
volScalarField D
(
 IOobject
 (
  "Darcy:D",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("0",K.dimensions(),0.0)
 );

volScalarField kc
(
 IOobject
 (
  "AbsorptionC",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedScalar("0",KC.dimensions(),0.0)
 );

volVectorField Uf
(
 IOobject
 (
  "Uf",
  mesh.time().timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 mesh,
 dimensionedVector("0",U.dimensions(),vector(0.0,0.0,0.0))
 );

//create the permeability 
labelList cellZoneIDs;
cellZoneIDs = mesh.cellZones().findIndices("Porous");
bool foundZone  = !cellZoneIDs.empty();
reduce(foundZone, orOp<bool>());
if(!foundZone && Pstream::master())
  {
    std::cout<< "cannot find porous cellZone Porous"<<std::endl;
  }

forAll(cellZoneIDs, zoneI)
{
  const labelList& cells = mesh.cellZones()[cellZoneIDs[zoneI]];
  forAll(cells, i)
    {
      const label cellI = cells[i];
      vector coords = mesh.C()[cellI];
      scalar a = coords[0]-xc.value();
      scalar b = coords[1]-yc.value();
      scalar r = Foam::sqrt(a*a+b*b);
      scalar ur = usrc.value()/r;
      D[cellI]  = K.value();
      kc[cellI] = KC.value();
      Uf[cellI] = vector(-ur*a/r, -ur*b/r, 0.0);
    }
}
