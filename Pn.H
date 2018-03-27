/*---------------------------------------------------------------------------*\
Class
    Foam::nongrayRadiation::Pn

Description
	Pn solver 

SourceFiles
    Pn.C

\*---------------------------------------------------------------------------*/

#ifndef nongrayRad_Pn_H
#define nongrayRad_Pn_H

#include "grayPn.H"
#include "basicNongrayRadiation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nongrayRadiation
{

/*---------------------------------------------------------------------------*\
                           Class Pn Declaration
\*---------------------------------------------------------------------------*/

class Pn
:
    public basicNongrayRadiation
{
    // Private data
        //- Configuration coefficients
        dictionary PnCoeffs_;
	
        //- spectral Pn
        PtrList<grayPn> Pns_;

        //- Radiative heat source
        autoPtr<volScalarField> Sh_;

        //- Incident radiation  [W/m2]
        PtrList<volScalarField> G_;

        //- Total radiative heat flux [W/m2]
        //volScalarField Qr_;

    // Private Member Functions
	//- Initialise
	void initialise();
	
        //- Disallow default bitwise copy construct
        Pn(const Pn&);

        //- Disallow default bitwise assignment
        void operator=(const Pn&);

public:

    //- Runtime type information
    TypeName("Pn");


    // Constructors

    //- Construct from components
    Pn
    (
        const volScalarField& T,
        const volScalarField& p,
        const basicMultiComponentMixture & C 
    );

    Pn
    (
        const thermoType thermo 
    );

    //- Construct from components
    Pn
    (
        const volScalarField& T,
        const volScalarField& p,
        const basicMultiComponentMixture & C ,
        const volScalarField & fv
    );

    Pn
    (
        const thermoType thermo ,
        const volScalarField & fv
    );


    //- Destructor
    virtual ~Pn();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

        // Return sources
		// source term to enthalpy
		virtual tmp<volScalarField> Sh() const;

		// source term to sensible enthalpy
		virtual tmp<volScalarField> Shs() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
