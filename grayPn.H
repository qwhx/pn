/*---------------------------------------------------------------------------*\
Class
    Foam::nongrayRadiation::grayPn

Description

    Gray Pn RTE solver. 
	
SourceFiles
    grayPn.C

\*---------------------------------------------------------------------------*/

#ifndef nongrayRad_grayPn_H
#define nongrayRad_grayPn_H

#include "PnFwds.H"
#include "coupleControl.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace nongrayRadiation
{

/*---------------------------------------------------------------------------*\
                           Class grayPn Declaration
\*---------------------------------------------------------------------------*/

class grayPn
{

    typedef coupleControl::solverPerformance solverPerformance;

    typedef coupleControl::solverPerformanceList solverPerformanceList;
    
    // Private data

        //- prefix
        const word prefix_;

        //- mesh
        const fvMesh & mesh_;

        //- Total absorption coefficient [1/m]
        const volScalarField & absc_;

        //- Total absorption coefficient [1/m]
        const volScalarField & scat_;

        //- Total emission [unit varies with spectral model]
        const volScalarField & E_;

        //- Configurations
        const dictionary & coeffs_;

        //- Order model
        autoPtr<PnOrderModel> order_;

        //- Write control
        Switch Write_;

        //- Couple control
        coupleControl couple_;

        //- List of pointers to radiative intensity
        PtrList<PnIntensity> Inm_;

        //- Incident radiation  [W/m2]
        volScalarField G_;

        //- Total radiative heat flux [W/m2]
        volScalarField Qr_;

        //- Emmited radiative heat flux [W/m2]
        volScalarField Qem_;

        //- Incidet radiative heat flux [W/m2]
        volScalarField Qin_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- Maximum omega weight
        scalar omegaMax_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        grayPn(const grayPn&);

        //- Disallow default bitwise assignment
        void operator=(const grayFvDOM&);

public:

    //- Runtime type information
    TypeName("grayPn");


    // Constructors

        //- Construct from components
        grayPn
        (
            const word prefix,
            const volScalarField& absc,
            const volScalarField& scat,
            const volScalarField& emis,
            const dictionary& dict
        );

    //- Destructor
    virtual ~grayPn();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Update G and calculate total heat flux on boundary
            void updateG();

        // Access

            //- Mesh
            inline const fvMesh& mesh() const;

            //- Intensity for order n, m
            inline const PnIntensity& Inm(const label n, const label m) const;

            //- Intensity for global order gid
            inline const PnIntensity& Inm(const label gid) const;

            //- Number of intensities
            inline label nInm() const;

            //- Const access to total absorption coefficient
            inline const volScalarField& absc() const;

            //- Const access to total scattering coefficient
            inline const volScalarField& scat() const;

            //- Const access to total emission
            inline const volScalarField& E() const;

            //- Const access to incident radiation field
            inline const volScalarField& G() const;

            //- Const access to total radiative heat flux field
            inline const volScalarField& Qr() const;

            //- Const access to incident radiative heat flux field
            inline const volScalarField& Qin() const;

            //- Const access to emitted radiative heat flux field
            inline const volScalarField& Qem() const;

            //- Const access to spherical harmonic coefficient model
            inline const PnOrderModel & order() const;

            //- Return omegaMax
            inline scalar omegaMax() const;

private:
            //- Check convergence
            bool checkConvergence(const solverPerformanceList &);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "grayPnI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
