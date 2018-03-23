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
	  // data for self calcuation
        //- Incident radiation / [W/m^2]
        volScalarField G_;

        //- Total Emission contribution [W/m^2]
        volScalarField E_;

        //- Radiative heat source [W/m^3]
        volScalarField Sh_;

    		//- Radiative heat flux to wall [W/m^2]
		    volScalarField Qr_;

		    //- internal calculation Gg_
		    PtrList<volScalarField> Gg_;

		    //- Extra coefficients
		    dictionary PnCoeffs_;

  	  	//- linear-anisotropic coefficient
		    scalar A1_; 

    // Private Member Functions

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
