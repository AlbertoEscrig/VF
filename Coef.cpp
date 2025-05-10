// =================================================================================================
// ======================================================================================== Coef.cpp
//
// Coeficientes del sistema
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Coef;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Tensor;
import :Malla;

import std;

namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// =========================================================================================== TCoef

template<std::size_t d, std::size_t r>
struct TCoef
{
  using TaN = std::array<double, TCelda<d>::MaxNCara>;

// ------------------------------------------------------------------------------------------- Datos

  double aP;
  TaN aN;
  TTensor<d, r> b;

// --------------------------------------------------------------------------------------- Funciones

private:
  template<std::size_t... i>
  TCoef
  constexpr Suma(std::index_sequence<i...>, TCoef const &Coef) const
    { return {aP + Coef.aP, {(aN[i] + Coef.aN[i])...}, b + Coef.b}; }

  template<std::size_t... i>
  TCoef
  constexpr Resta(std::index_sequence<i...>, TCoef const &Coef) const
    { return {aP - Coef.aP, {(aN[i] - Coef.aN[i])...}, b - Coef.b}; }

  template<std::size_t... i>
  TCoef
  constexpr Negativo(std::index_sequence<i...>) const
    { return {-aP, {(-aN[i])...}, -b}; }

  template<std::size_t... i>
  TCoef
  constexpr Producto(std::index_sequence<i...>, double const Escalar) const
    { return {aP * Escalar, {(aN[i] * Escalar)...}, b * Escalar}; }

// -------------------------------------------------------------------------------------- Operadores

public:
  TCoef
  constexpr operator +(TCoef const &Coef) const
    { return Suma(std::make_index_sequence<TCelda<d>::MaxNCara>{}, Coef); }

  TCoef
  constexpr operator -(TCoef const &Coef) const
    { return Resta(std::make_index_sequence<TCelda<d>::MaxNCara>{}, Coef); }

  TCoef
  constexpr operator -() const
    { return Negativo(std::make_index_sequence<TCelda<d>::MaxNCara>{}); }

  TCoef
  constexpr operator *(double const Escalar) const
    { return Producto(std::make_index_sequence<TCelda<d>::MaxNCara>{}, Escalar); }

  TCoef
  constexpr operator /(double const Escalar) const
    { return Producto(std::make_index_sequence<TCelda<d>::MaxNCara>{}, 1.0 / Escalar); }

// ------------------------------------------------------------------------------------------ Amigos

  TCoef
  friend constexpr operator *(double const Escalar, TCoef const &Coef)
    { return Coef * Escalar; }
};

} // VF
