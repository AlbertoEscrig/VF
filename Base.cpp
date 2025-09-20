// =================================================================================================
// ======================================================================================== Base.cpp
//
// Utilidades básicas
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Base;

// ========================================================================================= IMPORTS
// =================================================================================================

import std;

namespace VF
{
// ==================================================================== DECLARACIONES POR ADELANTADO
// ========================================================================================== TCelda

export
template<std::size_t>
class TCelda;

// =================================================================================================
// ========================================================================================== TMalla

export
template<std::size_t>
class TMalla;

// =================================================================================================
// ======================================================================================= TExprBase

template<typename>
class TExprBase;

// =================================================================================================
// ===================================================================================== TExprUnaria

template<std::size_t, std::size_t, typename, typename>
class TExprUnaria;

// =================================================================================================
// ==================================================================================== TExprBinaria

template<std::size_t, std::size_t, typename, typename, typename>
class TExprBinaria;

// =================================================================================================
// ========================================================================================== TCampo

export
template<std::size_t, std::size_t>
class TCampo;

// =================================================================================================
// ============================================================================================== Td

template<std::size_t, std::size_t>
class Td;

// =================================================================================================
// ======================================================================================== TSistema

export
template<std::size_t, std::size_t>
class TSistema;

// ================================================================= CONCEPTOS Y CONSTANTES EN LÍNEA
// =========================================================================================== CExpr

template<typename T>
concept CExpr = std::derived_from<T, TExprBase<T>>;

// =================================================================================================
// ========================================================================================= DimExpr

template<typename>
inline constexpr std::size_t DimExpr = std::numeric_limits<std::size_t>::max();

template<std::size_t d, std::size_t r, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
  requires CExpr<T<d, r, TRest...>>
inline constexpr std::size_t DimExpr<T<d, r, TRest...>> = d;

// =================================================================================================
// ======================================================================================= RangoExpr

template<typename>
inline constexpr std::size_t RangoExpr = std::numeric_limits<std::size_t>::max();

template<std::size_t d, std::size_t r, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
  requires CExpr<T<d, r, TRest...>>
inline constexpr std::size_t RangoExpr<T<d, r, TRest...>> = r;

// =================================================================================================
// ======================================================================================== CDimExpr

template<typename T, std::size_t d>
concept CDimExpr = CExpr<T> && DimExpr<T> == d;

// =================================================================================================
// ===================================================================================== CDimRanExpr

template<typename T, std::size_t d, std::size_t r>
concept CDimRanExpr = CDimExpr<T, d> && RangoExpr<T> == r;

// =================================================================================================
// ========================================================================================= EsCampo

template<typename>
inline constexpr bool EsCampo = false;

template<std::size_t d, std::size_t r>
inline constexpr bool EsCampo<TCampo<d, r>> = true;

// =================================================================================================
// ======================================================================================== EsTransi

template<typename>
inline constexpr bool EsTransi = false;

template<std::size_t d, std::size_t r>
inline constexpr bool EsTransi<Td<d, r>> = true;

template<std::size_t d, std::size_t r, typename T, typename TOp>
inline constexpr bool EsTransi<TExprUnaria<d, r, T, TOp>> = EsTransi<T>;

template<std::size_t d, std::size_t r, typename T, typename U, typename TOp>
inline constexpr bool EsTransi<TExprBinaria<d, r, T, U, TOp>> = EsTransi<T> || EsTransi<U>;

// ======================================================================================= FUNCIONES
// ============================================================================================= pow

export
template<std::unsigned_integral auto n, typename T>
decltype(auto)
constexpr pow(T const &x)
{
if constexpr (n == 0u)
  return static_cast<T>(1);
else if constexpr (n == 1u)
  return x;
else
  return x * pow<n - 1u>(x);
}

// =================================================================================================
// ============================================================================================ lerp

export
template<typename T, typename U>
auto
constexpr lerp(T const &lhs, U const &rhs, double const f)
  { return (1.0 - f) * lhs + f * rhs; }

} // VF
