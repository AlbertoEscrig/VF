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

constexpr std::size_t d = 2u;

constexpr VF::TTensor<d, 1u> i = {1.0, 0.0},
                             j = {0.0, 1.0};

constexpr VF::TTensor<d, 2u> I = i * i + j * j;

constexpr double H = 0.05;

constexpr double ν = 0.01,
                 Λ = 2.5,
                 ξ = 1.0,
                 C = 0.005;

constexpr VF::TTensor<d, 1u> U0 = 1.0 * i;

constexpr VF::TTensor<d, 2u> A0 = I / d;

constexpr double α1 = 0.6,
                 α2 = 0.3;

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

VF::TMalla<d>::Read("FolgarTucker.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampo<d, 0u> p;
VF::TCampo<d, 1u> U;
VF::TCampo<d, 2u> A;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TUInlet>("inlet");
A.DefCC<VF::TDirichlet>("inlet", A0);
U.DefCC<VF::TDirichlet>("wall");
p.DefCC<VF::TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0 * (1.0 - pow<2u>((VF::TMalla<d>::C() & j) / H));

A = A0;

while (true)
  {
// ------------------------------------------------------------------------------- Campos calculados

  auto const γ = (grad(U) + gradT(U)) / 2.0;
  auto const ω = (grad(U) - gradT(U)) / 2.0;

  VF::TCampo const τ = 2.0 * ν * Λ * (γ && A) * A;

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

  VF::TSistema const AEc =
    div(U * A) + 2.0 * ξ * (γ && A) * (+A) + C * d * mag(γ) * (+A)
    ==
    (A & ω) - (ω & A) + ξ * ((A & γ) + (γ & A)) + C * mag(γ) * I;

  solve(AEc, A, α2);
  }

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TMalla<d>::Write(ofs);
p.Write(ofs, "p");
U.Write(ofs, "U");
A.Write(ofs, "A");

return 0;
}
