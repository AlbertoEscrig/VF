// =================================================================================================
// ======================================================================================= Discr.cpp
//
// Discretización de operadores
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

//    EXPRESIÓN             OPERACIÓN         IMPLÍCITA
//    ---------             ---------         ---------
//    d(φ) / dt              ∂φ / ∂t              Sí
//
//      div(φ)                ∇ & φ               No
//
//    div(U * φ)           ∇ & (U * φ)            Sí
//
//     grad(φ)                ∇ * φ               No
//
//     gradT(φ)              (∇ * φ)ᵀ             No
//
//      lap(φ)              ∇ & ∇ * φ             Sí
//
// div(Γ * grad(φ))      ∇ & (Γ * ∇ * φ)          Sí
//
//  div(gradT(φ))          ∇ & (∇ * φ)ᵀ           No
//
//      rot(φ)                ∇ ^ φ               No
//
//      Sp(φ)                   φ                 Sí
//
//    Sp(ρ * φ)               ρ * φ               Sí

// ===================================================================== DECLARACIÓN DE LA PARTICIÓN
// =================================================================================================

export module VF:Discr;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;
import :Malla;
import :Coef;
import :Expr;
import :Campo;

import std;

namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// ============================================================================================== Td

template<std::size_t d, std::size_t r>
class Td : public TExprBase<Td<d, r>>
{
private:
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  Td(TCampo<d, r> const &φ_) :
    φ(φ_) {}

  template<bool>
  TCoef<d, r>
  Coef(TCelda<d> const &Celda) const
    { return {1.0, {}, -φ.Eval(Celda)}; }
};

// =================================================================================================
// ======================================================================================== TDivImpl

template<std::size_t d, std::size_t r, typename T>
class TDivImpl : public TExprBase<TDivImpl<d, r, T>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> U;
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

private:
  double
  SfUf(TCara<d> const &Cara) const requires CExpr<T>
    { return Cara.Sf & U.Eval(Cara); }

  double
  SfUf(TCara<d> const &Cara) const
    { return Cara.Sf & U; }

  template<bool>
  std::tuple<double, double, TTensor<d, r>>
  Coef(TCara<d> const &, double const) const;

public:
  TDivImpl(TExprBinaria<d, r + 1u, T, TCampo<d, r>, std::multiplies<>> const &Expr) :
    U(Expr.lhs), φ(Expr.rhs) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.Div(U, Celda); }

  template<bool>
  TCoef<d, r>
  Coef(TCelda<d> const &) const;
};

// =================================================================================================
// ======================================================================================== TDivExpl

template<std::size_t d, std::size_t r, CExpr T>
class TDivExpl : public TExprBase<TDivExpl<d, r, T>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TDivExpl(T const &Expr) :
    φ(Expr) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.Div(Celda); }

  template<bool>
  TTensor<d, r>
  Coef(TCelda<d> const &Celda) const
    { return Eval(Celda); }
};

// =================================================================================================
// =========================================================================================== TGrad

template<std::size_t d, std::size_t r, CExpr T>
class TGrad : public TExprBase<TGrad<d, r, T>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TGrad(T const &Expr) :
    φ(Expr) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.Grad(Celda); }

  template<bool>
  TTensor<d, r>
  Coef(TCelda<d> const &Celda) const
    { return Eval(Celda); }

// ------------------------------------------------------------------------------------------ Amigos

  friend class TLapBase<d, r - 1u>;
};

// =================================================================================================
// ========================================================================================== TGradT

template<std::size_t d, std::size_t r, CExpr T> requires (r == 2u)
class TGradT : public TExprBase<TGradT<d, r, T>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TGradT(T const &Expr) :
    φ(Expr) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.GradT(Celda); }

  template<bool>
  TTensor<d, r>
  Coef(TCelda<d> const &Celda) const
    { return Eval(Celda); }

// ------------------------------------------------------------------------------------------ Amigos

  friend class TLapT<d, r - 1u>;
};

// =================================================================================================
// ======================================================================================== TLapBase

template<std::size_t d, std::size_t r>
class TLapBase
{
protected:
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

protected:
  TLapBase(TCampo<d, r> const &φ_) :
    φ(φ_) {}

  TLapBase(TGrad<d, r + 1u, TCampo<d, r>> const &Expr) :
    φ(Expr.φ) {}

  std::tuple<double, double, TTensor<d, r>>
  Coef(TCara<d> const &, TTensor<d, r + 1u> const &) const;
};

// =================================================================================================
// ============================================================================================ TLap

template<std::size_t d, std::size_t r, typename T = void>
class TLap : public TLapBase<d, r>, public TExprBase<TLap<d, r, T>>
{
private:
  using TLapBase<d, r>::φ;
  std::conditional_t<EsCampo<T>, T const &, T const> Γ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TLap(TExprBinaria<d, r + 1u, T, TGrad<d, r + 1u, TCampo<d, r>>, std::multiplies<>> const &Expr) :
    TLapBase<d, r>(Expr.rhs), Γ(Expr.lhs) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.Lap(Γ, Celda); }

  template<bool>
  TCoef<d, r>
  Coef(TCelda<d> const &) const;
};

// =================================================================================================
// ================================================================================ TLap<d, r, void>

template<std::size_t d, std::size_t r>
class TLap<d, r, void> : public TLapBase<d, r>, public TExprBase<TLap<d, r, void>>
{
private:
  using TLapBase<d, r>::φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TLap(TCampo<d, r> const &φ_) :
    TLapBase<d, r>(φ_) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.Lap(Celda); }

  template<bool>
  TCoef<d, r>
  Coef(TCelda<d> const &) const;
};

// =================================================================================================
// =========================================================================================== TLapT

template<std::size_t d, std::size_t r>
class TLapT : public TExprBase<TLapT<d, r>>
{
private:
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TLapT(TGradT<d, r + 1u, TCampo<d, r>> const &Expr) :
    φ(Expr.φ) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.LapT(Celda); }

  template<bool>
  TTensor<d, r>
  Coef(TCelda<d> const &Celda) const
    { return Eval(Celda); }
};

// =================================================================================================
// ============================================================================================ TRot

template<std::size_t d, std::size_t r, CExpr T> requires (r == 1u)
class TRot : public TExprBase<TRot<d, r, T>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TRot(T const &Expr) :
    φ(Expr) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return φ.Rot(Celda); }

  template<bool>
  TTensor<d, r>
  Coef(TCelda<d> const &Celda) const
    { return Eval(Celda); }
};

// =================================================================================================
// ============================================================================================= TSp

template<std::size_t d, std::size_t r, typename T = void>
class TSp : public TExprBase<TSp<d, r, T>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> ρ;
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TSp(TExprBinaria<d, r, T, TCampo<d, r>, std::multiplies<>> const &Expr) :
    ρ(Expr.lhs), φ(Expr.rhs) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return ρ.Eval(Celda) * φ.Eval(Celda); }

  template<bool>
  TCoef<d, r>
  Coef(TCelda<d> const &) const;
};

// =================================================================================================
// ================================================================================= TSp<d, r, void>

template<std::size_t d, std::size_t r>
class TSp<d, r, void> : public TExprBase<TSp<d, r, void>>
{
private:
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TSp(TCampo<d, r> const &φ_) :
    φ(φ_) {}

  TTensor<d, r> const
  &Eval(TCelda<d> const &Celda) const
    { return φ.Eval(Celda); }

  template<bool>
  TCoef<d, r>
  Coef(TCelda<d> const &) const
    { return {1.0, {}, {}}; }
};

// ======================================================================== IMPLEMENTACIÓN DE CLASES
// ======================================================================================== TDivImpl

template<std::size_t d, std::size_t r, typename T>
template<bool EsTransi>
std::tuple<double, double, TTensor<d, r>>
TDivImpl<d, r, T>::Coef(TCara<d> const &Cara, double const SfUf) const
{
if constexpr (EsTransi)                          // Linear upwind
  {
  if (SfUf >= 0.0)
    return {SfUf, 0.0, SfUf * ((Cara.Cf - Cara.CeldaP().C) & φ.Grad(Cara.CeldaP()))};
  if (Cara.EsCC()) [[unlikely]]
    {
    auto const [aP, b] = φ.Coef(Cara);

    return {SfUf * aP, 0.0, SfUf * b};
    }
  return {0.0, SfUf, SfUf * ((Cara.Cf - Cara.CeldaN().C) & φ.Grad(Cara.CeldaN()))};
  }
else                                             // Upwind
  {
  if (SfUf >= 0.0)
    return {SfUf, 0.0, {}};
  if (Cara.EsCC()) [[unlikely]]
    {
    auto const [aP, b] = φ.Coef(Cara);

    return {SfUf * aP, 0.0, SfUf * b};
    }
  return {0.0, SfUf, {}};
  }
}

// =================================================================================================

template<std::size_t d, std::size_t r, typename T>
template<bool EsTransi>
TCoef<d, r>
TDivImpl<d, r, T>::Coef(TCelda<d> const &Celda) const
{
TCoef<d, r> Coef = {};

for (auto const &[i, Cara] : Celda | std::views::enumerate)
  {
  auto const [aP, aN, b] = TDivImpl::Coef<EsTransi>(Cara, SfUf(Cara));

  Coef.aP += aP;
  Coef.aN[i] = aN;
  Coef.b += b;
  }
return Coef / Celda.V;
}

// =================================================================================================
// ======================================================================================== TLapBase

template<std::size_t d, std::size_t r>
std::tuple<double, double, TTensor<d, r>>
TLapBase<d, r>::Coef(TCara<d> const &Cara, TTensor<d, r + 1u> const &gradφ) const
{
TVector<d> const &Sf = Cara.Sf;

if (Cara.EsCC()) [[unlikely]]
  {
  auto const [aP, b] = φ.GradCoef(Cara);

  return {mag(Sf) * aP, 0.0, mag(Sf) * b};
  }

TVector<d> const δ = Cara.δ(),
                 Sfδ = ((Sf & Sf) / (δ & Sf)) * δ;
double const aN = mag(Sfδ) / mag(δ);

return {-aN, aN, (Sf - Sfδ) & Cara.Interpola(gradφ, φ.Grad(Cara.CeldaN()))};
}

// =================================================================================================
// ============================================================================================ TLap

template<std::size_t d, std::size_t r, typename T>
template<bool>
TCoef<d, r>
TLap<d, r, T>::Coef(TCelda<d> const &Celda) const
{
TTensor<d, r + 1u> const gradφ = φ.Grad(Celda);
TCoef<d, r> Coef = {};

for (auto const &[i, Cara] : Celda | std::views::enumerate)
  {
  auto const [aP, aN, b] = TLapBase<d, r>::Coef(Cara, gradφ);
  double const Γf = Γ.Eval(Cara);

  Coef.aP += Γf * aP;
  Coef.aN[i] = Γf * aN;
  Coef.b += Γf * b;
  }
return Coef / Celda.V;
}

// =================================================================================================
// ================================================================================ TLap<d, r, void>

template<std::size_t d, std::size_t r>
template<bool>
TCoef<d, r>
TLap<d, r, void>::Coef(TCelda<d> const &Celda) const
{
TTensor<d, r + 1u> const gradφ = φ.Grad(Celda);
TCoef<d, r> Coef = {};

for (auto const &[i, Cara] : Celda | std::views::enumerate)
  {
  auto const [aP, aN, b] = TLapBase<d, r>::Coef(Cara, gradφ);

  Coef.aP += aP;
  Coef.aN[i] = aN;
  Coef.b += b;
  }
return Coef / Celda.V;
}

// =================================================================================================
// ============================================================================================= TSp

template<std::size_t d, std::size_t r, typename T>
template<bool>
TCoef<d, r>
TSp<d, r, T>::Coef(TCelda<d> const &Celda) const
{
if (double const ρP = ρ.Eval(Celda); ρP >= 0.0)
  return {ρP, {}, {}};
else
  return {0.0, {}, ρP * φ.Eval(Celda)};
}

// ======================================================================================= FUNCIONES
// =============================================================================================== d

export
template<std::size_t d1, std::size_t r>
Td<d1, r>
d(TCampo<d1, r> const &φ)
  { return {φ}; }

// =================================================================================================
// ============================================================================================= div

export
template<std::size_t d, std::size_t r, CDimRanExpr<d, 1u> T>
TDivImpl<d, r - 1u, T>
div(TExprBinaria<d, r, T, TCampo<d, r - 1u>, std::multiplies<>> const &Expr)
  { return {Expr}; }

// =================================================================================================

export
template<std::size_t d, std::size_t r, CDimRanExpr<d, 0u> T>
TLap<d, r - 1u, T>
div(TExprBinaria<d, r, T, TGrad<d, r, TCampo<d, r - 1u>>, std::multiplies<>> const &Expr)
  { return {Expr}; }

// =================================================================================================

export
template<std::size_t d>
TLapT<d, 1u>
div(TGradT<d, 2u, TCampo<d, 1u>> const &Expr)
  { return {Expr}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> >= 1u)
TDivExpl<DimExpr<T>, RangoExpr<T> - 1u, T>
div(T const &Expr)
  { return {Expr}; }

// =================================================================================================
// ============================================================================================ grad

export
template<CExpr T>
TGrad<DimExpr<T>, RangoExpr<T> + 1u, T>
grad(T const &Expr)
  { return {Expr}; }

// =================================================================================================
// =========================================================================================== gradT

export
template<CExpr T> requires (RangoExpr<T> == 1u)
TGradT<DimExpr<T>, 2u, T>
gradT(T const &Expr)
  { return {Expr}; }

// =================================================================================================
// ============================================================================================= lap

export
template<std::size_t d, std::size_t r>
TLap<d, r>
lap(TCampo<d, r> const &φ)
  { return {φ}; }

// =================================================================================================
// ============================================================================================= rot

export
template<CDimRanExpr<3u, 0u> T>
TRot<3u, 1u, T>
rot(T const &Expr)
  { return {Expr}; }

// =================================================================================================
// ============================================================================================== Sp

export
template<std::size_t d, std::size_t r>
TSp<d, r>
Sp(TCampo<d, r> const &φ)
  { return {φ}; }

// =================================================================================================

export
template<std::size_t d, std::size_t r, CDimRanExpr<d, 0u> T>
TSp<d, r, T>
Sp(TExprBinaria<d, r, T, TCampo<d, r>, std::multiplies<>> const &Expr)
  { return {Expr}; }

} // VF
