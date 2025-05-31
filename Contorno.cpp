// =================================================================================================
// ==================================================================================== Contorno.cpp
//
// Condiciones de contorno básicas
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Contorno;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;
import :Malla;
import :Coef;

import std;

export
namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// ================================================================================== TEsCCExplicita

template<typename>
struct TEsCCExplicita : std::false_type {};

// =================================================================================================
// ========================================================================================= TCCBase

template<std::size_t d, std::size_t r>
class TCCBase
{
public:
  virtual ~TCCBase() = default;

  std::tuple<double, TTensor<d, r>>
  virtual Coef(TCara<d> const &) const;

  std::tuple<double, TTensor<d, r>>
  virtual GradCoef(TCara<d> const &) const;
};

// =================================================================================================
// ====================================================================================== TDirichlet

template<std::size_t d, std::size_t r>
class TDirichlet : public TCCBase<d, r>
{
private:
  TTensor<d, r> const φb = {};                 // φf = φb

// --------------------------------------------------------------------------------------- Funciones

public:
  TDirichlet() = default;

  TDirichlet(TTensor<d, r> const &φb_) :
    φb(φb_) {}

  std::tuple<double, TTensor<d, r>>
  virtual Coef(TCara<d> const &) const override
    { return {0.0, φb}; }

  std::tuple<double, TTensor<d, r>>
  virtual GradCoef(TCara<d> const &Cara) const override
    { double const magL = mag(Cara.L()); return {-1.0 / magL, φb / magL}; }
};

// =================================================================================================
// ======================================================================================== TNeumann

template<std::size_t d, std::size_t r>
class TNeumann : public TCCBase<d, r>
{
private:
  TTensor<d, r> const gb = {};                // nf & ∇ * φf = gb

// --------------------------------------------------------------------------------------- Funciones

public:
  TNeumann() = default;

  TNeumann(TTensor<d, r> const &gb_) :
    gb(gb_) {}

  std::tuple<double, TTensor<d, r>>
  virtual Coef(TCara<d> const &Cara) const override
    { return {1.0, mag(Cara.L()) * gb}; }

  std::tuple<double, TTensor<d, r>>
  virtual GradCoef(TCara<d> const &) const override
    { return {0.0, gb}; }
};

// =================================================================================================
// ========================================================================================== TRobin

template<std::size_t d, std::size_t r>
class TRobin : public TCCBase<d, r>
{
private:
  double const k, h;
  TTensor<d, r> const φb;                     // k * (nf & ∇ * φf) + h * φf = h * φb

// --------------------------------------------------------------------------------------- Funciones

public:
  TRobin() = delete;

  TRobin(double const k_, double const h_, TTensor<d, r> const &φb_) :
    k(k_), h(h_), φb(φb_) {}

  std::tuple<double, TTensor<d, r>>
  virtual GradCoef(TCara<d> const &Cara) const override
    { double const magL = mag(Cara.L()); return {-1.0 / (k / h + magL), φb / (k / h + magL)}; }
};

// =================================================================================================
// ======================================================================================= TSimetria

template<std::size_t d, std::size_t r>
class TSimetria : public TCCBase<d, r>
{
private:
  TCampo<d, r> const &φ;

// --------------------------------------------------------------------------------------- Funciones

public:
  TSimetria(TCampo<d, r> const &φ_) :
    φ(φ_) {}

  std::tuple<double, TTensor<d, r>>
  virtual Coef(TCara<d> const &) const override;
};

template<std::size_t d, std::size_t r>
struct TEsCCExplicita<TSimetria<d, r>> : std::true_type {};

// =================================================================================================
// ================================================================================ TSimetria<d, 0u>

template<std::size_t d>
class TSimetria<d, 0u> : public TNeumann<d, 0u> {};

template<std::size_t d>
struct TEsCCExplicita<TSimetria<d, 0u>> : std::false_type {};

// =================================================================================================
// ====================================================================================== TPeriodica

template<std::size_t d, std::size_t r>
class TPeriodica : public TCCBase<d, r>
{
private:
  TCampo<d, r> const &φ;
  TTensor<d, 2u> const T;

// --------------------------------------------------------------------------------------- Funciones

public:
  TPeriodica() = delete;

  TPeriodica(TCampo<d, r> const &φ_, TTensor<d, 2u> const &T_ = TTensor<d, 2u>::I()) :
    φ(φ_), T(T_) {}

  std::tuple<double, TTensor<d, r>>
  virtual Coef(TCara<d> const &) const override;
};

template<std::size_t d, std::size_t r>
struct TEsCCExplicita<TPeriodica<d, r>> : std::true_type {};

// ======================================================================== IMPLEMENTACIÓN DE CLASES
// ========================================================================================= TCCBase

template<std::size_t d, std::size_t r>
std::tuple<double, TTensor<d, r>>
TCCBase<d, r>::Coef(TCara<d> const &Cara) const
{
auto const [aP, b] = GradCoef(Cara);
double const magL = mag(Cara.L());

return {1.0 + magL * aP, magL * b};
}

// =================================================================================================

template<std::size_t d, std::size_t r>
std::tuple<double, TTensor<d, r>>
TCCBase<d, r>::GradCoef(TCara<d> const &Cara) const
{
auto const [aP, b] = Coef(Cara);
double const magL = mag(Cara.L());

return {(aP - 1.0) / magL, b / magL};
}

// =================================================================================================
// ======================================================================================= TSimetria

template<std::size_t d, std::size_t r>
std::tuple<double, TTensor<d, r>>
TSimetria<d, r>::Coef(TCara<d> const &Cara) const
{
TVector<d> const &Sf = Cara.Sf;
TTensor<d, r> const &φP = φ.Eval(Cara.CeldaP());
TTensor<d, 2u> const T = TTensor<d, 2u>::I() - 2.0 * Sf * Sf / (Sf & Sf);

if constexpr (r == 1u)
  return {0.5, 0.5 * (T & φP)};
else if constexpr (r == 2u)
  return {0.5, 0.5 * (T & φP & T)};
else
  static_assert(false);
}

// =================================================================================================
// ====================================================================================== TPeriodica

template<std::size_t d, std::size_t r>
std::tuple<double, TTensor<d, r>>
TPeriodica<d, r>::Coef(TCara<d> const &Cara) const
{
TTensor<d, r> const &φN = φ.Eval(Cara.CeldaN());

if constexpr (r == 0u)
  return {0.5, 0.5 * φN};
else if constexpr (r == 1u)
  return {0.5, 0.5 * (T & φN)};
else if constexpr (r == 2u)
  return {0.5, 0.5 * (T & φN & T.T())};
else
  static_assert(false);
}

} // VF
