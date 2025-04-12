// =================================================================================================
// ====================================================================================== Tensor.cpp
//
// Tensor d-dimensional de rango r
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Tensor;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;

import std;

export
namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// =================================================================================== TTensor<d, r>

template<std::size_t d, std::size_t r1>
struct TTensor
{
  TTensor<d, r1 - 1u> TensorArr[d];

// --------------------------------------------------------------------------------------- Funciones

private:
  template<std::size_t... i>
  TTensor<d, r1>
  static constexpr Identidad(std::index_sequence<i...>)
    { return {(i / d == i % d ? 1.0 : 0.0)...}; }

  template<std::size_t... i>
  TTensor<d, r1>
  constexpr T(std::index_sequence<i...>) const
    { return {Slice<r1 - 1u, d, i>()...}; }

  template<std::size_t r2, std::size_t Stride, std::size_t Offset, std::size_t... i>
  TTensor<d, r2>
  constexpr Slice(std::index_sequence<i...>) const
    { return {std::bit_cast<std::array<double, pow<r1>(d)>>(*this)[i * Stride + Offset]...}; }

  template<std::size_t r2, std::size_t Stride, std::size_t Offset>
  TTensor<d, r2>
  constexpr Slice() const
    { return Slice<r2, Stride, Offset>(std::make_index_sequence<pow<r2>(d)>{}); }

  template<std::size_t... i>
  TTensor<d, r1>
  constexpr Suma(std::index_sequence<i...>, TTensor<d, r1> const &Tensor) const
    { return {(TensorArr[i] + Tensor[i])...}; }

  template<std::size_t... i>
  void
  constexpr SumaAsigna(std::index_sequence<i...>, TTensor<d, r1> const &Tensor)
    { ((TensorArr[i] += Tensor[i]), ...); }

  template<std::size_t... i>
  TTensor<d, r1>
  constexpr Resta(std::index_sequence<i...>, TTensor<d, r1> const &Tensor) const
    { return {(TensorArr[i] - Tensor[i])...}; }

  template<std::size_t... i>
  void
  constexpr RestaAsigna(std::index_sequence<i...>, TTensor<d, r1> const &Tensor)
    { ((TensorArr[i] -= Tensor[i]), ...); }

  template<std::size_t... i>
  TTensor<d, r1>
  constexpr Negativo(std::index_sequence<i...>) const
    { return {(-TensorArr[i])...}; }

  template<std::size_t r2, std::size_t... i>
  TTensor<d, r1 + r2>
  constexpr ProdExt(std::index_sequence<i...>, TTensor<d, r2> const &Tensor) const
    { return {(TensorArr[i] * Tensor)...}; }

  template<std::size_t... i>
  void
  constexpr ProdAsigna(std::index_sequence<i...>, double const Escalar)
    { ((TensorArr[i] *= Escalar), ...); }

  template<std::size_t... i>
  void
  constexpr DivAsigna(std::index_sequence<i...>, double const Escalar)
    { ((TensorArr[i] /= Escalar), ...); }

  template<std::size_t r2, std::size_t... i>
  TTensor<d, r1 + r2 - 2u>
  constexpr ProdInt(std::index_sequence<i...>, TTensor<d, r2> const &Tensor) const
    { return (... + (Slice<r1 - 1u, d, i>() * Tensor[i])); }

  template<std::size_t r2, std::size_t... i> requires (r1 == 1u)
  TTensor<d, r1 + r2 - 2u>
  constexpr ProdInt(std::index_sequence<i...>, TTensor<d, r2> const &Tensor) const
    { return (... + (TensorArr[i] * Tensor[i])); }

  template<std::size_t r2, std::size_t... i>
  TTensor<d, r1 + r2 - 4u>
  constexpr ProdIntDbl(std::index_sequence<i...>, TTensor<d, r2> const &Tensor) const
    { return (... + (Slice<r1 - 2u, d * d, i>() * Tensor.template Slice<r2 - 2u, 1u, i>())); }

  template<std::size_t... i> requires (r1 == 2u)
  TTensor<d, r1 - 2u>
  constexpr ProdIntDbl(std::index_sequence<i...>, TTensor<d, 2u> const &Tensor) const
    { return (... + (TensorArr[i] & Tensor[i])); }

  template<std::size_t... i>
  std::istream
  &ReadASCII(std::index_sequence<i...>, std::istream &is)
    { return (is >> ... >> TensorArr[i]); }

  template<std::size_t i, std::size_t... j> requires (r1 == 1u)
  std::ostream
  &WriteASCII(std::index_sequence<i, j...>, std::ostream &os) const
    { os << TensorArr[i]; return ((os << ' ' << TensorArr[j]), ...) << std::endl; }

  template<std::size_t... i>
  std::ostream
  &WriteASCII(std::index_sequence<i...>, std::ostream &os) const
    { return (os << ... << TensorArr[i]) << std::endl; }

public:
  TTensor<d, r1>
  static constexpr Identidad() requires (r1 == 2u)
    { return Identidad(std::make_index_sequence<d * d>{}); }

  TTensor<d, r1>
  constexpr T() const requires (r1 == 2u)
    { return T(std::make_index_sequence<d>{}); }

// -------------------------------------------------------------------------------------- Operadores

                                              // ------------------------------------ Suma (T1 + T2)
  TTensor<d, r1>
  constexpr operator +(TTensor<d, r1> const &Tensor) const
    { return Suma(std::make_index_sequence<d>{}, Tensor); }

                                              // ------------------------ Suma-asignación (T1 += T2)
  TTensor<d, r1>
  constexpr &operator +=(TTensor<d, r1> const &Tensor)
    { SumaAsigna(std::make_index_sequence<d>{}, Tensor); return *this; }

                                              // ----------------------------------- Resta (T1 - T2)
  TTensor<d, r1>
  constexpr operator -(TTensor<d, r1> const &Tensor) const
    { return Resta(std::make_index_sequence<d>{}, Tensor); }

                                              // ----------------------- Resta-asignación (T1 -= T2)
  TTensor<d, r1>
  constexpr &operator -=(TTensor<d, r1> const &Tensor)
    { RestaAsigna(std::make_index_sequence<d>{}, Tensor); return *this; }

                                              // ------------------------------------- Negativo (-T)
  TTensor<d, r1>
  constexpr operator -() const
    { return Negativo(std::make_index_sequence<d>{}); }

                                              // ----------------------- Producto exterior (T1 * T2)
  template<std::size_t r2>
  TTensor<d, r1 + r2>
  constexpr operator *(TTensor<d, r2> const &Tensor) const
    { return ProdExt(std::make_index_sequence<d>{}, Tensor); }

                                              // -------------- Producto-asignación escalar (T *= s)
  TTensor<d, r1>
  constexpr &operator *=(double const Escalar)
    { ProdAsigna(std::make_index_sequence<d>{}, Escalar); return *this; }

                                              // -------------- División-asignación escalar (T /= s)
  TTensor<d, r1>
  constexpr &operator /=(double const Escalar)
    { DivAsigna(std::make_index_sequence<d>{}, Escalar); return *this; }

                                              // ----------------------- Producto interior (T1 & T2)
  template<std::size_t r2>
  TTensor<d, r1 + r2 - 2u>
  constexpr operator &(TTensor<d, r2> const &Tensor) const
    { return ProdInt(std::make_index_sequence<d>{}, Tensor); }

                                              // ---------------- Producto interior doble (T1 && T2)
  template<std::size_t r2>
  TTensor<d, r1 + r2 - 4u>
  constexpr operator &&(TTensor<d, r2> const &Tensor) const
    { return ProdIntDbl(std::make_index_sequence<d * d>{}, Tensor); }

  TTensor<d, r1 - 2u>
  constexpr operator &&(TTensor<d, 2u> const &Tensor) const requires (r1 == 2u)
    { return ProdIntDbl(std::make_index_sequence<d>{}, Tensor); }

                                              // ------------ Acceso a componentes (T[i, j, k, ...])
  decltype(auto)
  constexpr operator [](this auto &Self, std::size_t const i)
    { return Self.TensorArr[i]; }

  decltype(auto)
  constexpr operator [](this auto &Self, std::size_t const i, std::integral auto const... j)
    { return Self[i][j...]; }

// ------------------------------------------------------------------------------------------ Amigos

  template<std::size_t, std::size_t>
  friend struct TTensor;

  std::istream
  friend &operator >>(std::istream &is, TTensor<d, r1> &Tensor)
    { return Tensor.ReadASCII(std::make_index_sequence<d>{}, is); }

  std::ostream
  friend &operator <<(std::ostream &os, TTensor<d, r1> const &Tensor)
    { return Tensor.WriteASCII(std::make_index_sequence<d>{}, os); }
};

// =================================================================================================
// ================================================================================== TTensor<d, 0u>

template<std::size_t d>
struct TTensor<d, 0u>
{
  double Data;

// --------------------------------------------------------------------------------------- Funciones

private:
  template<std::size_t r2, std::size_t... i>
  TTensor<d, r2>
  constexpr ProdExt(std::index_sequence<i...>, TTensor<d, r2> const &Tensor) const
    { return {(Data * Tensor[i])...}; }

public:
  constexpr TTensor() = default;

  constexpr TTensor(double const Escalar) :
    Data{Escalar} {}

// -------------------------------------------------------------------------------------- Operadores

  template<std::size_t r2>
  TTensor<d, r2>
  constexpr operator *(TTensor<d, r2> const &Tensor) const
    { return ProdExt(std::make_index_sequence<d>{}, Tensor); }

  TTensor<d, 0u>
  constexpr operator *(TTensor<d, 0u> const &Tensor) const
    { return Data * Tensor.Data; }

  constexpr operator double() const
    { return Data; }

  constexpr operator double &()
    { return Data; }
};

// =========================================================================================== ALIAS
// ====================================================================================== TVector<d>

template<std::size_t d>
using TVector = TTensor<d, 1u>;

// =================================================================================================
// ======================================================================================= TVector2D

using TVector2D = TVector<2u>;

// =================================================================================================
// ======================================================================================= TVector3D

using TVector3D = TVector<3u>;

// ================================================================= OPERADORES Y FUNCIONES GLOBALES
// =================================================================================== TTensor<d, r>

template<std::size_t r, std::size_t d>
TTensor<d, r>
constexpr operator *(TTensor<d, r> const &Tensor, double const Escalar)
  { return Tensor * TTensor<d, 0u>{Escalar}; }

// =================================================================================================

template<std::size_t r, std::size_t d>
TTensor<d, r>
constexpr operator *(double const Escalar, TTensor<d, r> const &Tensor)
  { return Tensor * TTensor<d, 0u>{Escalar}; }

// =================================================================================================

template<std::size_t r, std::size_t d>
TTensor<d, r>
constexpr operator /(TTensor<d, r> const &Tensor, double const Escalar)
  { return Tensor * TTensor<d, 0u>{1.0 / Escalar}; }

// =================================================================================================

template<std::size_t d>
double
constexpr mag(TTensor<d, 0u> const &Tensor)
  { return std::abs(Tensor); }

// =================================================================================================

template<std::size_t d>
double
constexpr mag(TTensor<d, 1u> const &Tensor)
  { return std::sqrt(Tensor & Tensor); }

// =================================================================================================

template<std::size_t d>
double
constexpr mag(TTensor<d, 2u> const &Tensor)
  { return std::sqrt(Tensor && Tensor); }

// =================================================================================================
// ======================================================================================= TVector3D

TVector3D
constexpr operator ^(TVector3D const &Vector1, TVector3D const &Vector2)
{
return {Vector1[1u] * Vector2[2u] - Vector1[2u] * Vector2[1u],
        Vector1[2u] * Vector2[0u] - Vector1[0u] * Vector2[2u],
        Vector1[0u] * Vector2[1u] - Vector1[1u] * Vector2[0u]};
}

} // VF
