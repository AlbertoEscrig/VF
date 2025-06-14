// =================================================================================================
// ======================================================================================== Expr.cpp
//
// Evaluación diferida de expresiones
//
// ================================================================= Copyright © 2025 Alberto Escrig
// ===================================================================== DECLARACIÓN DE LA PARTICIÓN

export module VF:Expr;

// ========================================================================================= IMPORTS
// =================================================================================================

import :Base;
import :Tensor;
import :Malla;
import :Coef;

import std;

namespace VF
{
// =========================================================================== DECLARACIÓN DE CLASES
// ======================================================================================= TExprBase

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
class TExprBase<T<d, r1, TRest...>>
{
private:
  using TSelf = T<d, r1, TRest...>;

// --------------------------------------------------------------------------------------- Funciones

public:
  auto
  Eval(TCelda<d> const &) const = delete;

  auto
  Eval(TCara<d> const &) const = delete;

  auto
  Eval(std::size_t const) const = delete;

  template<bool>
  auto
  Coef(TCelda<d> const &) const = delete;

  TTensor<d, r1 + 1u>
  Grad(this TSelf const &, TCelda<d> const &);

  TTensor<d, r1 + 1u>
  GradT(this TSelf const &, TCelda<d> const &) requires (r1 == 1u);

  TTensor<d, r1 - 1u>
  Div(this TSelf const &, TCelda<d> const &) requires (r1 >= 1u);

  TTensor<d, r1>
  Div(this TSelf const &, TTensor<d, 1u> const &, TCelda<d> const &);

  TTensor<d, r1>
  Div(this TSelf const &, CDimRanExpr<d, 1u> auto const &, TCelda<d> const &);

  TTensor<d, r1>
  Rot(this TSelf const &, TCelda<d> const &) requires (d == 3u && r1 == 1u);

// -------------------------------------------------------------------------------------- Operadores

  template<CDimRanExpr<d, r1> TExpr>
  TExprBinaria<d, r1, TSelf, TExpr, std::plus<>>
  operator +(this TSelf const &lhs, TExpr const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, r1, TSelf, TTensor<d, r1>, std::plus<>>
  operator +(this TSelf const &lhs, TTensor<d, r1> const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, 0u, TSelf, double, std::plus<>>
  operator +(this TSelf const &lhs, double const rhs) requires (r1 == 0u)
    { return {lhs, rhs}; }

  template<CDimRanExpr<d, r1> TExpr>
  TExprBinaria<d, r1, TSelf, TExpr, std::minus<>>
  operator -(this TSelf const &lhs, TExpr const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, r1, TSelf, TTensor<d, r1>, std::minus<>>
  operator -(this TSelf const &lhs, TTensor<d, r1> const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, 0u, TSelf, double, std::minus<>>
  operator -(this TSelf const &lhs, double const rhs) requires (r1 == 0u)
    { return {lhs, rhs}; }

  TExprUnaria<d, r1, TSelf, std::negate<>>
  operator -(this TSelf const &rhs)
    { return {rhs}; }

  template<CDimExpr<d> TExpr, std::size_t r2 = RangoExpr<TExpr>>
  TExprBinaria<d, r1 + r2, TSelf, TExpr, std::multiplies<>>
  operator *(this TSelf const &lhs, TExpr const &rhs)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r1 + r2, TSelf, TTensor<d, r2>, std::multiplies<>>
  operator *(this TSelf const &lhs, TTensor<d, r2> const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, r1, TSelf, double, std::multiplies<>>
  operator *(this TSelf const &lhs, double const rhs)
    { return {lhs, rhs}; }

  template<CDimRanExpr<d, 0u> TExpr>
  TExprBinaria<d, r1, TSelf, TExpr, std::divides<>>
  operator /(this TSelf const &lhs, TExpr const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, r1, TSelf, double, std::divides<>>
  operator /(this TSelf const &lhs, double const rhs)
    { return {lhs, rhs}; }

  template<CDimExpr<d> TExpr, std::size_t r2 = RangoExpr<TExpr>>
  TExprBinaria<d, r1 + r2 - 2u, TSelf, TExpr, std::bit_and<>>
  operator &(this TSelf const &lhs, TExpr const &rhs) requires (r1 >= 1u && r2 >= 1u)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r1 + r2 - 2u, TSelf, TTensor<d, r2>, std::bit_and<>>
  operator &(this TSelf const &lhs, TTensor<d, r2> const &rhs) requires (r1 >= 1u && r2 >= 1u)
    { return {lhs, rhs}; }

  template<CDimExpr<d> TExpr, std::size_t r2 = RangoExpr<TExpr>>
  TExprBinaria<d, r1 + r2 - 4u, TSelf, TExpr, std::logical_and<>>
  operator &&(this TSelf const &lhs, TExpr const &rhs) requires (r1 >= 2u && r2 >= 2u)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r1 + r2 - 4u, TSelf, TTensor<d, r2>, std::logical_and<>>
  operator &&(this TSelf const &lhs, TTensor<d, r2> const &rhs) requires (r1 >= 2u && r2 >= 2u)
    { return {lhs, rhs}; }

  template<CDimRanExpr<d, r1> TExpr>
  TExprBinaria<3u, 1u, TSelf, TExpr, std::bit_xor<>>
  operator ^(this TSelf const &lhs, TExpr const &rhs) requires (d == 3u && r1 == 1u)
    { return {lhs, rhs}; }

  TExprBinaria<3u, 1u, TSelf, TTensor<3u, 1u>, std::bit_xor<>>
  operator ^(this TSelf const &lhs, TTensor<3u, 1u> const &rhs) requires (d == 3u && r1 == 1u)
    { return {lhs, rhs}; }

  TSistema<d, r1>
  operator ==(this TSelf const &lhs, CDimRanExpr<d, r1> auto const &rhs)
    { return {lhs, rhs}; }

  TSistema<d, r1>
  operator ==(this TSelf const &lhs, TTensor<d, r1> const &rhs)
    { return {lhs, rhs}; }

  decltype(auto)
  operator [](this TSelf const &Self, std::size_t const i) requires requires { Self.Eval(i); }
    { return Self.Eval(i); }

  decltype(auto)
  operator [](this TSelf const &Self, std::size_t const i)
    { return Self.Eval(TMalla<d>::Celda(i)); }

// ------------------------------------------------------------------------------------------ Amigos

  TExprBinaria<d, r1, TTensor<d, r1>, TSelf, std::plus<>>
  friend operator +(TTensor<d, r1> const &lhs, TSelf const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, 0u, double, TSelf, std::plus<>>
  friend operator +(double const lhs, TSelf const &rhs) requires (r1 == 0u)
    { return {lhs, rhs}; }

  TExprBinaria<d, r1, TTensor<d, r1>, TSelf, std::minus<>>
  friend operator -(TTensor<d, r1> const &lhs, TSelf const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, 0u, double, TSelf, std::minus<>>
  friend operator -(double const lhs, TSelf const &rhs) requires (r1 == 0u)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r1 + r2, TTensor<d, r2>, TSelf, std::multiplies<>>
  friend operator *(TTensor<d, r2> const &lhs, TSelf const &rhs)
    { return {lhs, rhs}; }

  TExprBinaria<d, r1, double, TSelf, std::multiplies<>>
  friend operator *(double const lhs, TSelf const &rhs)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r2, TTensor<d, r2>, TSelf, std::divides<>>
  friend operator /(TTensor<d, r2> const &lhs, TSelf const &rhs) requires (r1 == 0u)
    { return {lhs, rhs}; }

  TExprBinaria<d, 0u, double, TSelf, std::divides<>>
  friend operator /(double const lhs, TSelf const &rhs) requires (r1 == 0u)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r1 + r2 - 2u, TTensor<d, r2>, TSelf, std::bit_and<>>
  friend operator &(TTensor<d, r2> const &lhs, TSelf const &rhs) requires (r1 >= 1u && r2 >= 1u)
    { return {lhs, rhs}; }

  template<std::size_t r2>
  TExprBinaria<d, r1 + r2 - 4u, TTensor<d, r2>, TSelf, std::logical_and<>>
  friend operator &&(TTensor<d, r2> const &lhs, TSelf const &rhs) requires (r1 >= 2u && r2 >= 2u)
    { return {lhs, rhs}; }

  TExprBinaria<3u, 1u, TTensor<3u, 1u>, TSelf, std::bit_xor<>>
  friend operator ^(TTensor<3u, 1u> const &lhs, TSelf const &rhs) requires (d == 3u && r1 == 1u)
    { return {lhs, rhs}; }
};

// =================================================================================================
// ===================================================================================== TExprUnaria

template<std::size_t d, std::size_t r, typename T, typename TOp>
class TExprUnaria : public TExprBase<TExprUnaria<d, r, T, TOp>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> rhs;
  [[no_unique_address]] TOp const Op = {};

// --------------------------------------------------------------------------------------- Funciones

public:
  TExprUnaria(T const &rhs_) :
    rhs(rhs_) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const
    { return Op(rhs.Eval(Celda)); }

  TTensor<d, r>
  Eval(TCara<d> const &Cara) const
    { return Op(rhs.Eval(Cara)); }

  TTensor<d, r>
  Eval(std::size_t const i) const
    { return Op(rhs[i]); }

  template<bool EsTransi>
  decltype(auto)
  Coef(TCelda<d> const &Celda) const
    { return Op(rhs.template Coef<EsTransi>(Celda)); }
};

// =================================================================================================
// ==================================================================================== TExprBinaria

template<std::size_t d, std::size_t r, typename T, typename U, typename TOp>
class TExprBinaria : public TExprBase<TExprBinaria<d, r, T, U, TOp>>
{
private:
  std::conditional_t<EsCampo<T>, T const &, T const> lhs;
  std::conditional_t<EsCampo<U>, U const &, U const> rhs;
  [[no_unique_address]] TOp const Op = {};

// --------------------------------------------------------------------------------------- Funciones

public:
  TExprBinaria(T const &lhs_, U const &rhs_) :
    lhs(lhs_), rhs(rhs_) {}

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const requires (CExpr<T> && CExpr<U>)
    { return Op(lhs.Eval(Celda), rhs.Eval(Celda)); }

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const requires CExpr<T>
    { return Op(lhs.Eval(Celda), rhs); }

  TTensor<d, r>
  Eval(TCelda<d> const &Celda) const requires CExpr<U>
    { return Op(lhs, rhs.Eval(Celda)); }

  TTensor<d, r>
  Eval(TCara<d> const &) const requires (CExpr<T> && CExpr<U>);

  TTensor<d, r>
  Eval(TCara<d> const &Cara) const requires CExpr<T>
    { return Op(lhs.Eval(Cara), rhs); }

  TTensor<d, r>
  Eval(TCara<d> const &Cara) const requires CExpr<U>
    { return Op(lhs, rhs.Eval(Cara)); }

  TTensor<d, r>
  Eval(std::size_t const i) const requires (CExpr<T> && CExpr<U>)
    { return Op(lhs[i], rhs[i]); }

  TTensor<d, r>
  Eval(std::size_t const i) const requires CExpr<T>
    { return Op(lhs[i], rhs); }

  TTensor<d, r>
  Eval(std::size_t const i) const requires CExpr<U>
    { return Op(lhs, rhs[i]); }

  template<bool EsTransi>
  decltype(auto)
  Coef(TCelda<d> const &Celda) const requires (CExpr<T> && CExpr<U>)
    { return Op(lhs.template Coef<EsTransi>(Celda), rhs.template Coef<EsTransi>(Celda)); }

  template<bool EsTransi>
  decltype(auto)
  Coef(TCelda<d> const &Celda) const requires CExpr<T>
    { return Op(lhs.template Coef<EsTransi>(Celda), rhs); }

  template<bool EsTransi>
  decltype(auto)
  Coef(TCelda<d> const &Celda) const requires CExpr<U>
    { return Op(lhs, rhs.template Coef<EsTransi>(Celda)); }

// ------------------------------------------------------------------------------------------ Amigos

  friend class TDivImpl<d, r - 1u, T>;

  friend class TLap<d, r - 1u, T>;

  friend class TSp<d, r, T>;
};

// ======================================================================== IMPLEMENTACIÓN DE CLASES
// ======================================================================================= TExprBase

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
TTensor<d, r1 + 1u>
TExprBase<T<d, r1, TRest...>>::Grad(this TSelf const &Self, TCelda<d> const &Celda)
{
TTensor<d, r1 + 1u> gradφ = {};

for (auto &Cara : Celda)
  gradφ += Cara.Sf * Self.Eval(Cara);
return gradφ / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
TTensor<d, r1 + 1u>
TExprBase<T<d, r1, TRest...>>::GradT(this TSelf const &Self,
                                     TCelda<d> const &Celda) requires (r1 == 1u)
{
TTensor<d, r1 + 1u> gradφT = {};

for (auto &Cara : Celda)
  gradφT += Self.Eval(Cara) * Cara.Sf;
return gradφT / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
TTensor<d, r1 - 1u>
TExprBase<T<d, r1, TRest...>>::Div(this TSelf const &Self,
                                   TCelda<d> const &Celda) requires (r1 >= 1u)
{
TTensor<d, r1 - 1u> divφ = {};

for (auto &Cara : Celda)
  divφ += Cara.Sf & Self.Eval(Cara);
return divφ / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
TTensor<d, r1>
TExprBase<T<d, r1, TRest...>>::Div(this TSelf const &Self, TTensor<d, 1u> const &U,
                                   TCelda<d> const &Celda)
{
TTensor<d, r1> divφ = {};

for (auto &Cara : Celda)
  divφ += (Cara.Sf & U) * Self.Eval(Cara);
return divφ / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
TTensor<d, r1>
TExprBase<T<d, r1, TRest...>>::Div(this TSelf const &Self, CDimRanExpr<d, 1u> auto const &U,
                                   TCelda<d> const &Celda)
{
TTensor<d, r1> divφ = {};

for (auto &Cara : Celda)
  divφ += (Cara.Sf & U.Eval(Cara)) * Self.Eval(Cara);
return divφ / Celda.V;
}

// =================================================================================================

template<std::size_t d, std::size_t r1, typename... TRest,
         template<std::size_t, std::size_t, typename...> typename T>
TTensor<d, r1>
TExprBase<T<d, r1, TRest...>>::Rot(this TSelf const &Self,
                                   TCelda<d> const &Celda) requires (d == 3u && r1 == 1u)
{
TTensor<d, r1> rotφ = {};

for (auto &Cara : Celda)
  rotφ += Cara.Sf ^ Self.Eval(Cara);
return rotφ / Celda.V;
}

// =================================================================================================
// ==================================================================================== TExprBinaria

template<std::size_t d, std::size_t r, typename T, typename U, typename TOp>
TTensor<d, r>
TExprBinaria<d, r, T, U, TOp>::Eval(TCara<d> const &Cara) const requires (CExpr<T> && CExpr<U>)
{
if (Cara.EsCC()) [[unlikely]]
  return Op(lhs.Eval(Cara), rhs.Eval(Cara));
return Cara.Interpola(Eval(Cara.CeldaP()), Eval(Cara.CeldaN()));
}

} // VF
