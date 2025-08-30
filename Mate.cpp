// =================================================================================================
// ======================================================================================== Mate.cpp
//
// Funciones matemáticas que se aplican a expresiones
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Mate;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;
import :Malla;
import :Expr;

import std;

namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// ============================================================================================ TMin

struct TMin
{
  double
  static inline operator ()(double const lhs, double const rhs)
    { return std::min(lhs, rhs); }
};

// =================================================================================================
// ============================================================================================ TMax

struct TMax
{
  double
  static inline operator ()(double const lhs, double const rhs)
    { return std::max(lhs, rhs); }
};

// =================================================================================================
// ============================================================================================ TMag

template<std::size_t d, std::size_t r>
struct TMag
{
  double
  static inline operator ()(TTensor<d, r> const &Tensor)
    { return mag(Tensor); }
};

// ======================================================================================= FUNCIONES
// ============================================================================================= min

export
template<CExpr T, CExpr U>
  requires (DimExpr<T> == DimExpr<U> && RangoExpr<T> == 0u && RangoExpr<U> == 0u)
TExprBinaria<DimExpr<T>, 0u, T, U, TMin>
min(T const &lhs, U const &rhs)
  { return {lhs, rhs}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> == 0u)
TExprBinaria<DimExpr<T>, 0u, double, T, TMin>
min(double const lhs, T const &rhs)
  { return {lhs, rhs}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> == 0u)
TExprBinaria<DimExpr<T>, 0u, T, double, TMin>
min(T const &lhs, double const rhs)
  { return {lhs, rhs}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> == 0u)
double
min(T const &Expr)
{
double Min = std::numeric_limits<double>::max();

#pragma omp parallel for reduction(min : Min)
for (std::size_t i = 0u; i < TMalla<DimExpr<T>>::NCelda(); ++i)
  if (double const Val = Expr[i]; Val < Min)
    Min = Val;
return Min;
}

// =================================================================================================
// ============================================================================================= max

export
template<CExpr T, CExpr U>
  requires (DimExpr<T> == DimExpr<U> && RangoExpr<T> == 0u && RangoExpr<U> == 0u)
TExprBinaria<DimExpr<T>, 0u, T, U, TMax>
max(T const &lhs, U const &rhs)
  { return {lhs, rhs}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> == 0u)
TExprBinaria<DimExpr<T>, 0u, double, T, TMax>
max(double const lhs, T const &rhs)
  { return {lhs, rhs}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> == 0u)
TExprBinaria<DimExpr<T>, 0u, T, double, TMax>
max(T const &lhs, double const rhs)
  { return {lhs, rhs}; }

// =================================================================================================

export
template<CExpr T> requires (RangoExpr<T> == 0u)
double
max(T const &Expr)
{
double Max = std::numeric_limits<double>::min();

#pragma omp parallel for reduction(max : Max)
for (std::size_t i = 0u; i < TMalla<DimExpr<T>>::NCelda(); ++i)
  if (double const Val = Expr[i]; Val > Max)
    Max = Val;
return Max;
}

// =================================================================================================
// ============================================================================================= sum

export
template<CExpr T, std::size_t d = DimExpr<T>, std::size_t r = RangoExpr<T>>
TTensor<d, r>
sum(T const &Expr)
{
TTensor<d, r> Sum = {};

#pragma omp declare reduction(+ : TTensor<d, r> : omp_out += omp_in) initializer(omp_priv = {})
#pragma omp parallel for reduction(+ : Sum)
for (std::size_t i = 0u; i < TMalla<d>::NCelda(); ++i)
  Sum += Expr[i];
return Sum;
}

// =================================================================================================
// ============================================================================================= mag

export
template<CExpr T, std::size_t d = DimExpr<T>, std::size_t r = RangoExpr<T>>
TExprUnaria<d, 0u, T, TMag<d, r>>
mag(T const &Expr)
  { return {Expr}; }

// =================================================================================================
// ================================================================================= FN_MATE_ESCALAR

#define FN_MATE_ESCALAR(FN)                                                                        \
  struct T##FN                                                                                     \
  {                                                                                                \
    double                                                                                         \
    static inline operator ()(double const Escalar)                                                \
      { return std::FN(Escalar); }                                                                 \
  };                                                                                               \
                                                                                                   \
  export                                                                                           \
  template<CExpr T, std::size_t d = DimExpr<T>> requires (RangoExpr<T> == 0u)                      \
  TExprUnaria<d, 0u, T, T##FN>                                                                     \
  FN(T const &Expr)                                                                                \
    { return {Expr}; }

FN_MATE_ESCALAR(sqrt)
FN_MATE_ESCALAR(cbrt)
FN_MATE_ESCALAR(exp)
FN_MATE_ESCALAR(log)
FN_MATE_ESCALAR(sin)
FN_MATE_ESCALAR(cos)
FN_MATE_ESCALAR(tan)
FN_MATE_ESCALAR(asin)
FN_MATE_ESCALAR(acos)
FN_MATE_ESCALAR(atan)

} // VF
