#include "grayPn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace nongrayRadiation
    {
        defineTypeNameAndDebug(grayPn, 0);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nongrayRadiation::grayPn::initialise()
{
    order_ = PnOrderModel::New(coeffs_);
    label nInm = order_->size();
    Inm_.setSize(nInm);
    forAll (Inm_, i)
    {
        label n = order_->n(i);
        label m = order_->m(i);
        word prefix (prefix_ +"I"+ Foam::name(i));
        Inm_.set
        (
            i,
            new PnIntensity
            (
                i,
                prefix,
                *this,
                absc_,
                scat_,
                E_,
                coeffs_
            )
        );
    }
    Info<< "grayPn : Allocated " << Inm_.size() 
        << " intensity fields." << nl;

    forAll (Inm__, i)
    {
        Inm_[i].updateCoupledIntensityBoundaryConditions();
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nongrayRadiation::grayPn::grayPn
(
    const word prefix,
    const volScalarField& absc,
    const volScalarField& scat,
    const volScalarField& emis,
    const dictionary& dict
)
:
    prefix_(prefix),
    mesh_(absc.mesh()),
    absc_(absc),
    scat_(scat),
    E_(emis),
    coeffs_(dict),
    Write_(coeffs_.lookup("Write")),
    order_(),
    couple_(absc.mesh(), "Pn"),
    Inm_(),
    G_
    (
        IOobject
        (
            prefix_+":G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            bool(Write_) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    Qr_
    (
        IOobject
        (
            prefix_+":Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            bool(Write_) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    Qem_
    (
        IOobject
        (
            prefix_+":Qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            bool(Write_) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qem", dimMass/pow3(dimTime), 0.0)
    ),
    Qin_
    (
        IOobject
        (
            prefix_+":Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            bool(Write_) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimMass/pow3(dimTime), 0.0)
    ),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    omegaMax_(0)
{
    initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nongrayRadiation::grayPn::~grayPn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nongrayRadiation::grayPn::read()
{
        // Only reading solution parameters - not changing intensity orders
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
}


void Foam::nongrayRadiation::grayPn::calculate()
{
    // Update all intensities
    scalar maxResidual = 0.0;
    label radIter = 0;
    bool nextLoop = true;
    while (couple_.loop())
    {
        Info<< "Radiation solver iter: " << radIter << endl;
        radIter++;
        maxResidual = 0.0;
        solverPerformanceList sol(nInm());
        forAll(Inm_, intensityI)
        {
            sol[intensityI] = Inm_[intensityI].correct();
        }
        nextLoop = checkConvergence(sol);
    } 

    updateG();
}


void Foam::nongrayRadiation::grayPn::updateG()
{
    G_ = Inm_[order_->gid(0,0)]*4*pi;
}

bool Foam::nongrayRadiation::grayPn::checkConvergence
(
    const grayPn::solverPerformanceList & sol
)
{
    scalar maxInitRes = 0.0;
    scalar maxFinalRes = 0.0;
    label maxInitResIntensityId = -1;
    label maxFinalResIntensityId = -1;

    forAll (sol, intensityI)
    {
        const scalar intensityInitRes = sol[intensityI].initialResidual();
        const scalar intensityFinalRes = sol[intensityI].finalResidual();
        if (intensityInitRes > maxInitRes)
        {
            maxInitRes = intensityInitRes;
            maxInitResIntensityId = intensityI;
        }
        if (intensityFinalRes > maxFinalRes)
        {
            maxFinalRes = intensityFinalRes;
            maxFinalResIntensityId = intensityI;
        }
    }

    Info<< "Intensity "<< maxInitResIntensityId<< " has max initial residual: "<< "\t"<< maxInitRes << endl;
    Info<< "Intensity "<< maxFinalResIntensityId<<" has max final residual: "<< "\t"<< maxFinalRes << endl;

    couple_.checkConvergence(sol);
    return !couple_.toExit();
}
// ************************************************************************* //
