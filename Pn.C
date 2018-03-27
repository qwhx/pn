#include "Pn.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace nongrayRadiation
    {
        defineTypeNameAndDebug(Pn, 0);
        addToRunTimeSelectionTable
        (
            basicNongrayRadiation,
            Pn,
            composition
        );

        addToRunTimeSelectionTable
        (
            basicNongrayRadiation,
            Pn,
            thermo
        );

        addToRunTimeSelectionTable
        (
            basicNongrayRadiation,
            Pn,
            compositionSoot
        );

        addToRunTimeSelectionTable
        (
            basicNongrayRadiation,
            Pn,
            thermoSoot
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nongrayRadiation::Pn::initialise()
{
    Pns_.resize(spectralField().size());
    // Construct absorption field for each wavelength
    forAll(Pns_, PnI)
    {
		    word prefix("Pn::Quad"+Foam::name(PnI));
		    const basicSpectralField & sp = spectralPtr_();
		    const volScalarField & absc = sp.k(PnI);
		    const volScalarField & scat = sp.scat(PnI);
		    const volScalarField & emis = sp.E(PnI);

        Pns_.set
        (
            PnI,
            new grayPn
            (
                prefix,
				        absc,
				        scat,
				        emis,
				        PnCoeffs_
            )
        );
    }

	  const word ShName("Pn::Sh");
	  Sh_.set
	  (
		    new volScalarField
		    (
			      IOobject
			      (
				        ShName,
				        mesh_.time().timeName(),
				        mesh_,
				        IOobject::NO_READ,
        		    bool(Write_) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
			      ),
			      mesh_,
            dimensionedScalar("Sh", dimMass/dimLength/pow3(dimTime), 0.0)
		    )
	  );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
Foam::nongrayRadiation::Pn::Pn
(
	const volScalarField& T,
	const volScalarField& p,
	const basicMultiComponentMixture & C
)
:
    basicNongrayRadiation(T,p,C),
	  PnCoeffs_(this->subDict("PnCoeffs")),
	  Pns_(),
	  Sh_(),
	  G_()
{
    initialise();
}

Foam::nongrayRadiation::Pn::Pn
(
    const basicNongrayRadiation::thermoType thermo
)
:
    basicNongrayRadiation(thermo),
    PnCoeffs_(this->subDict("PnCoeffs")),
    Pns_(),
    Sh_(),
    G_()
{
    initialise();
}

Foam::nongrayRadiation::Pn::Pn
(
    const volScalarField& T,
    const volScalarField& p,
    const basicMultiComponentMixture & C,
    const volScalarField & fv
)
:
    basicNongrayRadiation(T,p,C, fv),
    PnCoeffs_(this->subDict("PnCoeffs")),
    Pns_(),
    Sh_(),
    G_()
{
    initialise();
}

Foam::nongrayRadiation::Pn::Pn
(
    const basicNongrayRadiation::thermoType thermo,
    const volScalarField & fv
)
:
    basicNongrayRadiation(thermo, fv),
    PnCoeffs_(this->subDict("PnCoeffs")),
    Pns_(),
    Sh_(),
    G_()
{
    initialise();
}
} 
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nongrayRadiation::Pn::~Pn()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nongrayRadiation::Pn::read()
{
    // Reserved for future use
        
    if (basicNongrayRadiation::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::nongrayRadiation::Pn::calculate()
{
    Sh_()=dimensionedScalar("zero",Sh_->dimensions(), 0.0);

    spectralPtr_->correct();

    Info<<"Foam::nongrayRadiation::Pn::calculate() : "
		    <<"Finished update spectral variables. "
		    <<"ClockTime = " << Sh_->mesh().time().elapsedClockTime() 
		    << " s" << nl << endl;

    forAll (Pns_, PnI)
    {
        Pns_[PnI].calculate();
        Pns_[PnI].updateG();
        
        const scalar w = spectralPtr_->W(PnI);
        const volScalarField & Gi = Pns_[PnI].G();
        const volScalarField & Ei = Pns_[PnI].E();
        const volScalarField & ki = Pns_[PnI].absc();
        Sh_() += w*ki*(Gi - Ei);
    }
}

Foam::tmp<Foam::volScalarField> Foam::nongrayRadiation::Pns::Sh () const
{
    return Foam::tmp<volScalarField>(Sh_());
}

Foam::tmp<Foam::volScalarField> Foam::nongrayRadiation::Pns::Shs () const
{
    return Foam::tmp<volScalarField>(Sh_());
}
// ************************************************************************* //
