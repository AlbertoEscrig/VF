// =================================================================================================
// ======================================================================================= Sistema.h
//
// Resolución del sistema por el método del gradiente biconjugado estabilizado (BiCGSTAB)
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Sistema;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;
import :Malla;
import :Coef;
import :Expr;
import :Mate;
import :Campo;

import std;

namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// ======================================================================================== TSistema

export
template<std::size_t d, std::size_t r>
class TSistema
{
private:
  template<std::size_t = d, std::size_t = r>
  class TΣaN : public TExprBase<TΣaN<d, r>>
  {
  private:
    TSistema const &Sistema;
    TCampo<d, r> const &φ;

  public:
    TΣaN(TSistema const &Sistema_, TCampo<d, r> const &φ_) :
      Sistema(Sistema_), φ(φ_) {}

    TTensor<d, r>
    Eval(TCelda<d> const &Celda) const
      { return Sistema.ΣaN(Celda.ID, φ); }

    TTensor<d, r>
    Eval(std::size_t const i) const
      { return Sistema.ΣaN(i, φ); }
  };

// ------------------------------------------------------------------------------------------- Datos

private:
  static constexpr double ε = 1e-8;             // Error relativo a ||b||
  static constexpr std::size_t MaxIt = 10'000u; // Número máximo de iteraciones

  std::vector<typename TCoef<d, r>::TaN> aNVec; // Términos no diagonales

public:
  TCampo<d, 0u> aP;                             // Diagonal
  TCampo<d, r> b;                               // Términos independientes

// --------------------------------------------------------------------------------------- Funciones

private:
  double
  static Dot(TCampo<d, 0u> const &lhs, TCampo<d, 0u> const &rhs)
    { return sum(lhs * rhs); }

  double
  static Dot(TCampo<d, 1u> const &lhs, TCampo<d, 1u> const &rhs)
    { return sum(lhs & rhs); }

  double
  static Dot(TCampo<d, 2u> const &lhs, TCampo<d, 2u> const &rhs)
    { return sum(lhs && rhs); }

  TTensor<d, r>
  ΣaN(std::size_t const, TCampo<d, r> const &) const;

  void
  Solve(TCampo<d, r> &, double const) const;

public:
  TSistema() = delete;

  template<CExpr T, typename U>
  TSistema(T const &, U const &);

  void
  DefRef(std::size_t const ID, TTensor<d, r> const &Ref)
    { b[ID] += aP[ID] * Ref; aP[ID] *= 2.0; }

  TΣaN<>
  ΣaN(TCampo<d, r> const &φ) const
    { return {*this, φ}; }

// ------------------------------------------------------------------------------------------ Amigos

  friend class TΣaN<>;

  void
  friend solve(TSistema const &Sistema, TCampo<d, r> &Campo, double const f = 1.0)
    { Sistema.Solve(Campo, f); }
};

// ======================================================================== IMPLEMENTACIÓN DE CLASES
// ======================================================================================== TSistema

template<std::size_t d, std::size_t r>
template<CExpr T, typename U>
TSistema<d, r>::TSistema(T const &lhs, U const &rhs) :
  aNVec(TMalla<d>::NCelda())
{
#pragma omp parallel for
for (std::size_t i = 0u; i < TMalla<d>::NCelda(); ++i)
  {
  TCelda<d> const &Celda = TMalla<d>::Celda(i);
  TCoef const Coef = lhs.template Coef<EsTransi<T>>(Celda);

  aP[i] = Coef.aP;
  aNVec[i] = Coef.aN;
  if constexpr (CExpr<U>)
    b[i] = rhs.Eval(Celda) - Coef.b;
  else
    b[i] = rhs - Coef.b;
  }
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TSistema<d, r>::ΣaN(std::size_t const i, TCampo<d, r> const &φ) const
{
TTensor<d, r> aNφ = {};

for (auto const &[aN, ID] : std::views::zip(aNVec[i], TMalla<d>::IDVec(i)))
  if (ID != MaxID) [[likely]]
    aNφ += aN * φ[ID];
return aNφ;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
void
TSistema<d, r>::Solve(TCampo<d, r> &φ, double const f) const
{
double const tol2 = pow<2u>(ε) * Dot(b, b);
TCampo<d, r> x = φ,
             e = b - aP * x - ΣaN(x),
             e0 = e,
             p, s,
             y, z,
             v, t;
double ρ = 1.0,
       α = 1.0,
       ω = 1.0,
       e0e0 = Dot(e0, e0);

for (std::size_t It = 0u; It < MaxIt; ++It)
  {
  double const ρOld = ρ;

  if (Dot(e, e) < tol2)
    break;
  ρ = Dot(e0, e);
  if (std::abs(ρ) < pow<2u>(std::numeric_limits<double>::epsilon()) * e0e0)
    {
    e0 = e = b + aP * (1.0 - f) / f * φ - aP / f * x - ΣaN(x);
    ρ = e0e0 = Dot(e0, e0);
    }
  p = e + (ρ / ρOld) * (α / ω) * (p - ω * v);
  y = p / aP;
  v = aP * y / f + ΣaN(y);
  α = ρ / Dot(e0, v);
  s = e - α * v;
  z = s / aP;
  t = aP * z / f + ΣaN(z);
  ω = Dot(t, s) / Dot(t, t);
  x += α * y + ω * z;
  e = s - ω * t;
  }
φ = std::move(x);
}

} // VF
