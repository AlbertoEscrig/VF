// =================================================================================================
// ================================================================================ FolgarTucker.cpp
//
// Orientación de una suspesión de fibras que circula por un canal
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

// -------------------------------------------------------------------------------------- Constantes

constexpr VF::TTensor<2u, 1u> i = {1.0, 0.0},
                              j = {0.0, 1.0};

constexpr VF::TTensor<2u, 2u> I = i * i + j * j;

constexpr double H = 0.05;

constexpr double ν = 0.01,
                 Λ = 2.5,
                 ξ = 1.0,
                 C = 0.005;

constexpr VF::TTensor<2u, 1u> U0 = 1.0 * i;

constexpr VF::TTensor<2u, 2u> A0 = 0.5 * I;

constexpr double α1 = 0.5,
                 α2 = 0.2;

// ------------------------------------------------------ Condición de contorno para U en la entrada

template<std::size_t d, std::size_t r> requires (r == 1u)
class TUInlet : public VF::TCCBase<d, r>
{
public:
  std::tuple<double, VF::TTensor<d, r>>
  virtual Coef(VF::TCara<d> const &Cara) const override
    { return {0.0, U0 * (1.0 - VF::pow<2u>((Cara.Cf & j) / H))}; }
};

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla<2u>::Read("FolgarTucker.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampo<2u, 0u> p;
VF::TCampo<2u, 1u> U;
VF::TCampo<2u, 2u> A;
VF::TCampo<2u, 2u> τ;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TUInlet>("inlet");
A.DefCC<VF::TDirichlet>("inlet", A0);
U.DefCC<VF::TDirichlet>("wall");
τ.DefCC<VF::TDirichlet>("wall");
p.DefCC<VF::TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0 * (1.0 - pow<2u>((VF::TMalla<2u>::C() & j) / H));

A = A0;

while (true)
  {
// ---------------------------------------------------------------------------------- Momento lineal

  VF::TSistema const UEc = div(U * U) - ν * lap(U) == div(τ) - grad(p);

  solve(UEc, U, α1);

// ------------------------------------------------------------------------------------- Continuidad

  U = VF::TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

  VF::TCampo const pOld = p;

  solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 1e-4 * sum(mag(p)))
    break;

  p = lerp(p, pOld, α1);

// --------------------------------------------------------------------------- Tensor de orientación

  auto const γ = (grad(U) + gradT(U)) / 2.0;
  auto const ω = (grad(U) - gradT(U)) / 2.0;

  VF::TSistema const AEc =
    div(U * A) + 2.0 * ξ * (γ && A) * (+A) + 3.0 * C * mag(γ) * (+A)
    ==
    (A & ω) - (ω & A) + ξ * ((γ & A) + (A & γ)) + C * mag(γ) * I;

  solve(AEc, A, α2);

  τ = 2.0 * ν * Λ * (γ && A) * A;
  }

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TMalla<2u>::Write(ofs);
p.Write(ofs, "p");
U.Write(ofs, "U");
A.Write(ofs, "A");

return 0;
}
