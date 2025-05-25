// =================================================================================================
// ======================================================================================= Campo.cpp
//
// Campo tensorial
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Campo;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;
import :Malla;
import :Coef;
import :Contorno;
import :Expr;

import std;

export
namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// ========================================================================================== TCampo

template<std::size_t d, std::size_t r>
class TCampo : public TExprBase<TCampo<d, r>>
{
private:
  using TCCPtr = std::unique_ptr<TCCBase<d, r>>;

// ------------------------------------------------------------------------------------------- Datos

private:
  std::size_t const   NCelda     = TMalla<d>::NCelda();
  TTensor<d, r>       *TensorPtr = new TTensor<d, r>[NCelda]{};
  std::vector<TCCPtr> CCPtrVec   = CCNeumann() | std::views::take(TMalla<d>::NCC())
                                               | std::ranges::to<std::vector>();

// --------------------------------------------------------------------------------------- Funciones

private:
  std::generator<TCCPtr>
  static CCNeumann()                             // Por defecto Neumann homogénea
    { while (true) co_yield std::make_unique<TNeumann<d, r>>(); }

  template<typename T, typename... TArgs>
  void
  DefCC(std::false_type, std::string_view const CCStr, TArgs &&...Args)
    { CCPtrVec[TMalla<d>::CCID(CCStr)] = std::make_unique<T>(std::forward<TArgs>(Args)...); }

  template<typename T, typename... TArgs>
  void
  DefCC(std::true_type, std::string_view const CCStr, TArgs &&...Args)
    { DefCC<T>(std::false_type{}, CCStr, *this, std::forward<TArgs>(Args)...); }

  TTensor<d, r>
  Lap(TCelda<d> const &, TCara<d> const &, TTensor<d, r + 1u> const &) const;

  TTensor<d, r>
  LapT(TCelda<d> const &, TCara<d> const &, TTensor<d, r + 1u> const &) const;

public:
  TCampo() = default;

  ~TCampo()
    { delete[] TensorPtr; }

  TCampo(TCampo const &Campo)
    { Asigna(Campo); }

  TCampo(TCampo &&Campo) :
    TensorPtr(std::exchange(Campo.TensorPtr, nullptr)) {}

  TCampo(TTensor<d, r> const &Tensor)
    { Asigna(Tensor); }

  TCampo(CDimRanExpr<d, r> auto const &Expr)
    { Asigna(Expr); }

  void
  Asigna(TCampo const &Campo)
    { std::copy_n(Campo.TensorPtr, NCelda, TensorPtr); }

  void
  Asigna(TTensor<d, r> const &Tensor)
    { std::fill_n(TensorPtr, NCelda, Tensor); }

  void
  Asigna(CDimRanExpr<d, r> auto const &);

  void
  Asigna(std::string_view const, TTensor<d, r> const &Tensor);

  void
  Asigna(std::string_view const, CDimRanExpr<d, r> auto const &);

  template<template<std::size_t, std::size_t> typename T, typename... TArgs>
    requires std::derived_from<T<d, r>, TCCBase<d, r>>
  void
  DefCC(std::string_view const CCStr, TArgs &&...Args)
    { DefCC<T<d, r>>(TEsCCExplicita<T<d, r>>{}, CCStr, std::forward<TArgs>(Args)...); }

  TTensor<d, r> const
  &Eval(TCelda<d> const &Celda) const
    { return TensorPtr[Celda.ID]; }

  TTensor<d, r>
  Eval(TCara<d> const &) const;

  TTensor<d, r> const
  &Eval(std::size_t const i) const
    { return TensorPtr[i]; }

  template<bool>
  TTensor<d, r> const
  &Coef(TCelda<d> const &Celda) const
    { return Eval(Celda); }

  std::tuple<double, TTensor<d, r>>
  Coef(TCara<d> const &Cara) const
    { return CCPtrVec[Cara.CCID]->Coef(Cara); }

  using TExprBase<TCampo<d, r>>::Grad;

  std::tuple<double, TTensor<d, r>>
  GradCoef(TCara<d> const &Cara) const
    { return CCPtrVec[Cara.CCID]->GradCoef(Cara); }

  using TExprBase<TCampo<d, r>>::GradT;

  TTensor<d, r>
  Lap(TCelda<d> const &) const;

  TTensor<d, r>
  Lap(CDimRanExpr<d, 0u> auto const &, TCelda<d> const &) const;

  TTensor<d, r>
  LapT(TCelda<d> const &) const requires (r == 1u);

  void
  Write(std::ostream &, std::string_view const) const;

// -------------------------------------------------------------------------------------- Operadores

  TCampo
  &operator =(TCampo const &Campo)
    { Asigna(Campo); return *this; }

  TCampo
  &operator =(TCampo &&Campo)
    { delete[] TensorPtr; TensorPtr = std::exchange(Campo.TensorPtr, nullptr); return *this; }

  TCampo
  &operator =(TTensor<d, r> const &Tensor)
    { Asigna(Tensor); return *this; }

  TCampo
  &operator =(CDimRanExpr<d, r> auto const &Expr)
    { Asigna(Expr); return *this; }

  TCampo
  &operator +=(TTensor<d, r> const &);

  TCampo
  &operator +=(CDimRanExpr<d, r> auto const &);

  TCampo
  &operator -=(TTensor<d, r> const &);

  TCampo
  &operator -=(CDimRanExpr<d, r> auto const &);

  using TExprBase<TCampo<d, r>>::operator [];

  TTensor<d, r>
  &operator [](std::size_t const i)
    { return TensorPtr[i]; }

// ------------------------------------------------------------------------------------------ Amigos

  TTensor<d, r> const
  friend *begin(TCampo const &Campo)
    { return Campo.TensorPtr; }

  TTensor<d, r>
  friend *begin(TCampo &Campo)
    { return Campo.TensorPtr; }

  TTensor<d, r> const
  friend *end(TCampo const &Campo)
    { return Campo.TensorPtr + Campo.NCelda; }

  TTensor<d, r>
  friend *end(TCampo &Campo)
    { return Campo.TensorPtr + Campo.NCelda; }
};

// ===================================================================== GUÍA DE DEDUCCIÓN EXPLÍCITA
// ========================================================================================== TCampo

template<CExpr T>
TCampo(T const &) -> TCampo<DimExpr<T>, RangoExpr<T>>;

// =========================================================================================== ALIAS
// =================================================================================== TCampoEscalar

template<std::size_t d>
using TCampoEscalar = TCampo<d, 0u>;

// =================================================================================================
// ================================================================================= TCampoEscalar2D

using TCampoEscalar2D = TCampoEscalar<2u>;

// =================================================================================================
// ================================================================================= TCampoEscalar3D

using TCampoEscalar3D = TCampoEscalar<3u>;

// =================================================================================================
// ================================================================================= TCampoVectorial

template<std::size_t d>
using TCampoVectorial = TCampo<d, 1u>;

// =================================================================================================
// =============================================================================== TCampoVectorial2D

using TCampoVectorial2D = TCampoVectorial<2u>;

// =================================================================================================
// =============================================================================== TCampoVectorial3D

using TCampoVectorial3D = TCampoVectorial<3u>;

// ======================================================================== IMPLEMENTACIÓN DE CLASES
// ========================================================================================== TCampo

template<std::size_t d, std::size_t r>
void
TCampo<d, r>::Asigna(CDimRanExpr<d, r> auto const &Expr)
{
#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  TensorPtr[i] = Expr[i];
}

// =================================================================================================

template<std::size_t d, std::size_t r>
void
TCampo<d, r>::Asigna(std::string_view const GrpStr, TTensor<d, r> const &Tensor)
{
std::size_t const GrpID = TMalla<d>::GrpID(GrpStr);

#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  if (TMalla<d>::Celda(i).GrpID == GrpID)
    TensorPtr[i] = Tensor;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
void
TCampo<d, r>::Asigna(std::string_view const GrpStr, CDimRanExpr<d, r> auto const &Expr)
{
std::size_t const GrpID = TMalla<d>::GrpID(GrpStr);

#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  if (TMalla<d>::Celda(i).GrpID == GrpID)
    TensorPtr[i] = Expr[i];
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TCampo<d, r>
&TCampo<d, r>::operator +=(TTensor<d, r> const &Tensor)
{
#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  TensorPtr[i] += Tensor;
return *this;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TCampo<d, r>
&TCampo<d, r>::operator +=(CDimRanExpr<d, r> auto const &Expr)
{
#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  TensorPtr[i] += Expr[i];
return *this;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TCampo<d, r>
&TCampo<d, r>::operator -=(TTensor<d, r> const &Tensor)
{
#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  TensorPtr[i] -= Tensor;
return *this;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TCampo<d, r>
&TCampo<d, r>::operator -=(CDimRanExpr<d, r> auto const &Expr)
{
#pragma omp parallel for
for (std::size_t i = 0u; i < NCelda; ++i)
  TensorPtr[i] -= Expr[i];
return *this;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TCampo<d, r>::Eval(TCara<d> const &Cara) const
{
if (Cara.EsCC()) [[unlikely]]
  {
  auto const [aP, b] = Coef(Cara);

  return aP * Eval(Cara.CeldaP()) + b;
  }
return Cara.Interpola(Eval(Cara.CeldaP()), Eval(Cara.CeldaN()));
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TCampo<d, r>::Lap(TCelda<d> const &Celda, TCara<d> const &Cara,
                  TTensor<d, r + 1u> const &gradφ) const
{
if (Cara.EsCC()) [[unlikely]]
  {
  auto const [aP, b] = GradCoef(Cara);

  return mag(Cara.Sf) * (aP * Eval(Celda) + b);
  }
return Cara.Sf & Cara.Interpola(gradφ, Grad(Cara.CeldaN()));
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TCampo<d, r>::Lap(TCelda<d> const &Celda) const
{
TTensor<d, r + 1u> const gradφ = Grad(Celda);
TTensor<d, r> lapφ = {};

for (auto &Cara : Celda)
  lapφ += Lap(Celda, Cara, gradφ);
return lapφ / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TCampo<d, r>::Lap(CDimRanExpr<d, 0u> auto const &Γ, TCelda<d> const &Celda) const
{
TTensor<d, r + 1u> const gradφ = Grad(Celda);
TTensor<d, r> lapφ = {};

for (auto &Cara : Celda)
  lapφ += Γ.Eval(Cara) * Lap(Celda, Cara, gradφ);
return lapφ / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TCampo<d, r>::LapT(TCelda<d> const &Celda, TCara<d> const &Cara,
                   TTensor<d, r + 1u> const &gradφT) const
{
if (Cara.EsCC()) [[unlikely]]
  {
  TVector<d> const nf = Cara.nf();
  auto const [aP, b] = GradCoef(Cara);

  return Cara.Sf & (gradφT + (aP * Eval(Celda) + b - (gradφT & nf)) * nf);
  }
return Cara.Sf & Cara.Interpola(gradφT, GradT(Cara.CeldaN()));
}

// =================================================================================================

template<std::size_t d, std::size_t r>
TTensor<d, r>
TCampo<d, r>::LapT(TCelda<d> const &Celda) const requires (r == 1u)
{
TTensor<d, r + 1u> const gradφT = GradT(Celda);
TTensor<d, r> lapφT = {};

for (auto &Cara : Celda)
  lapφT += LapT(Celda, Cara, gradφT);
return lapφT / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r>
void
TCampo<d, r>::Write(std::ostream &os, std::string_view const Nombre) const
{
if constexpr (r == 0u)
  {
  os << "SCALARS " << Nombre << " double 1" << std::endl;
  os << "LOOKUP_TABLE default" << std::endl;
  for (auto &Tensor : *this)
    os << Tensor << std::endl;
  }
else if constexpr (r == 1u)
  {
  os << "VECTORS " << Nombre << " double" << std::endl;
  for (auto &Tensor : *this)
    if constexpr (d == 2u)
      os << Tensor[0u] << ' ' << Tensor[1u] << " 0" << std::endl;
    else if constexpr (d == 3u)
      os << Tensor;
  }
else if constexpr (r == 2u)
  {
  os << "TENSORS " << Nombre << " double" << std::endl;
  for (auto &Tensor : *this)
    if constexpr (d == 2u)
      {
      os << Tensor[0u, 0u] << ' ' << Tensor[0u, 1u] << " 0" << std::endl;
      os << Tensor[1u, 0u] << ' ' << Tensor[1u, 1u] << " 0" << std::endl;
      os << "0 0 0" << std::endl;
      }
    else if constexpr (d == 3u)
      os << Tensor;
  }
}

} // VF
