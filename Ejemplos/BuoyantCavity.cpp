// =================================================================================================
// =============================================================================== BuoyantCavity.cpp
//
// Cavidad con corrientes de flotación
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

// ------------------------------------------------------ Condición de contorno para p - ρ * (g & C)

template<std::size_t d, std::size_t r> requires (r == 0u)
class TGradPresionDin : public VF::TCCBase<d, r>
{
private:
  VF::TCampoEscalar<d> const &ρ;
  VF::TVector<d> const g;

public:
  TGradPresionDin(VF::TCampoEscalar<d> const &ρ_, VF::TVector<d> const &g_) :
    ρ(ρ_), g(g_) {}

  std::tuple<double, VF::TTensor<d, r>>
  virtual GradCoef(VF::TCara<d> const &Cara) const override
  {
  auto const [aP, b] = ρ.GradCoef(Cara);

  return {0.0, -(g & Cara.Cf) * (aP * ρ.Eval(Cara.CeldaP()) + b)};
  }
};

// -------------------------------------------------------------------------------------- Constantes

constexpr VF::TVector2D g = {0.0, -9.81};

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

VF::TMalla2D::Read("BuoyantCavity.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoVectorial2D U;
VF::TCampoEscalar2D p,
                    ρ,
                    T;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<VF::TDirichlet>("hot");
p.DefCC<TGradPresionDin>("hot", ρ, g);
T.DefCC<VF::TDirichlet>("hot", Thot);
ρ.DefCC<VF::TDirichlet>("hot", p0 * M / (R * Thot));

U.DefCC<VF::TDirichlet>("cold");
p.DefCC<TGradPresionDin>("cold", ρ, g);
T.DefCC<VF::TDirichlet>("cold", Tcold);
ρ.DefCC<VF::TDirichlet>("cold", p0 * M / (R * Tcold));

U.DefCC<VF::TDirichlet>("adiabatic");
p.DefCC<TGradPresionDin>("adiabatic", ρ, g);

// --------------------------------------------------------------------------- Condiciones iniciales

T = T0;

// ------------------------------------------------------------------------------- Campos calculados

auto const C = VF::TMalla2D::C();

while (true)
  {
  ρ = p0 * M / (R * T);

// ---------------------------------------------------------------------------------- Momento lineal

  VF::TSistema const UEc =
    div(ρ * U * U) - μ * lap(U) == μ * grad(div(U)) / 3.0 - grad(p) - (g & C) * grad(ρ);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  VF::TCampo const φ = ρ * (UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP; // CC Neumann homogéneas

  VF::TCampo const pOld = p;

  VF::TSistema pEc = div((ρ / UEc.aP) * grad(p)) == div(φ);

  pEc.DefRef(0u, 0.0);

  solve(pEc, p);

  U = φ / ρ - grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 1e-4)
    break;

  p = lerp(p, pOld, α);

// ----------------------------------------------------------------------------------------- Energía

  solve(div(ρ * U * T) - div(ρ * U) * (+T) - μ / Pr * lap(T) == 0, T);
  }

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");
T.Write(ofs, "T");
ρ.Write(ofs, "rho");

return 0;
}
