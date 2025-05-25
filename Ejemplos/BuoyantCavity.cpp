// =================================================================================================
// =============================================================================== BuoyantCavity.cpp
//
// Cavidad con corrientes de flotación
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// ------------------------------------------------------ Condición de contorno para p - ρ * (g & C)

template <std::size_t d, std::size_t r> requires (r == 0u)
class TGradPresionDin : public TCCBase<d, r>
{
private:
  TCampoEscalar<d> const &ρ;
  TVector<d> const g;

public:
  TGradPresionDin() = delete;

  TGradPresionDin(TCampoEscalar<d> const &ρ_, TVector<d> const &g_) :
    ρ(ρ_), g(g_) {}

  std::tuple<double, TTensor<d, r>>
  virtual GradCoef(TCara<d> const &) const override;
};

template <std::size_t d, std::size_t r> requires (r == 0u)
std::tuple<double, TTensor<d, r>>
TGradPresionDin<d, r>::GradCoef(TCara<d> const &Cara) const
{
auto const [aP, b] = ρ.GradCoef(Cara);

return {0.0, -(g & Cara.Cf) * (aP * ρ.Eval(Cara.CeldaP()) + b)};
}

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector2D g = {0.0, -9.81};

constexpr double μ  = 1.831e-5, // Pa s
                 M  = 28.96,    // kg / kmol
                 R  = 8314.5,   // J / kmol K
                 Pr = 0.705;

constexpr double p0 = 1e5,
                 T0 = 293.0,
                 Thot  = T0 + 2.0,
                 Tcold = T0 - 2.0;

constexpr double α = 0.7;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("BuoyantCavity.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p,
                ρ,
                T;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("hotWall");
p.DefCC<TGradPresionDin>("hotWall", ρ, g);
T.DefCC<TDirichlet>("hotWall", Thot);
ρ.DefCC<TDirichlet>("hotWall", p0 * M / (R * Thot));

U.DefCC<TDirichlet>("coldWall");
p.DefCC<TGradPresionDin>("coldWall", ρ, g);
T.DefCC<TDirichlet>("coldWall", Tcold);
ρ.DefCC<TDirichlet>("coldWall", p0 * M / (R * Tcold));

U.DefCC<TDirichlet>("adiabaticWall");
p.DefCC<TGradPresionDin>("adiabaticWall", ρ, g);

// --------------------------------------------------------------------------- Condiciones iniciales

T = T0;

// ------------------------------------------------------------------------------- Campos calculados

auto const C = TMalla2D::C();

while (true)
  {
  ρ = p0 * M / (R * T);

// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc =
    div(ρ * U * U) - μ * lap(U) == μ * div(gradT(U)) / 3.0 - grad(p) - (g & C) * grad(ρ);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  TCampo const φ = ρ * (UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP; // CC Neumann homogéneas

  TCampo const pOld = p;

  TSistema pEc = div((ρ / UEc.aP) * grad(p)) == div(φ);

  pEc.DefRef(0u, 0.0);

  solve(pEc, p);

  U = φ / ρ - grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 1e-4)
    break;

  p = lerp(p, pOld, α);

// ----------------------------------------------------------------------------------------- Energía

  solve(div(ρ * U * T) - div(ρ * U) * Sp(T) - μ / Pr * lap(T) == 0, T);
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");
T.Write(ofs, "T");
ρ.Write(ofs, "rho");

return 0;
}
