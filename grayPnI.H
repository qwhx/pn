inline const Foam::fvMesh &
Foam::nongrayRadiation::grayPn::mesh() const
{
    return G().mesh();
}

inline const Foam::nongrayRadiation::PnIntensity&
Foam::nongrayRadiation::grayPn::Inm(const label gid) const
{
    return Inm_[gid];
}

inline const Foam::nongrayRadiation::PnIntensity&
Foam::nongrayRadiation::grayPn::Inm(const label n, const label m) const
{
    return Inm_[order_->gid(n,m)];
}

inline Foam::label Foam::nongrayRadiation::grayPn::nInm() const
{
    return Inm_.size();
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::absc() const
{
    return absc_;
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::scat() const
{
    return scat_;
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::E() const
{
    return E_;
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::G() const
{
    return G_;
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::Qr() const
{
    return Qr_;
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::Qin() const
{
    return Qin_;
}

inline const Foam::volScalarField& 
Foam::nongrayRadiation::grayPn::Qem() const
{
    return Qem_;
}

inline const Foam::nongrayRadiation::directionalQuadratureModel &
Foam::nongrayRadiation::grayPn::order() const
{
    return order_();
}

inline Foam::scalar Foam::nongrayRadiation::grayPn::omegaMax() const
{
    return omegaMax_;
}

// ************************************************************************* //
