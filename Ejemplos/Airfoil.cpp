// =================================================================================================
// ===================================================================================== Airfoil.cpp
//
// Flujo en torno a un perfil aerodinámico
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double Cμ = 0.09,
                 κ  = 0.41,
                 β  = 3.0 / 40.0,
                 γ  = 5.0 / 9.0,
                 σk = 0.5,
                 σω = 0.5;

constexpr VF::TVector2D U0 = {10.0, 0.0};

constexpr double k0 = 0.015,
                 ω0 = 100.0;

constexpr double ν = 1.5e-5;

constexpr double α = 0.7;

constexpr std::size_t NCorNoOrto = 3u;

// -------------------------------------------------------------------- Condición de contorno para ω

template<std::size_t d, std::size_t r> requires (r == 0u)
class TCCω : public VF::TCCBase<d, 0u>
{
public:
  std::tuple<double, VF::TTensor<d, 0u>>
  virtual Coef(VF::TCara<d> const &Cara) const override
    { double const L = Cara.L(); return {0.0, 6.0 * ν / (β * L * L)}; }
};

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla2D::Read("NACA2412.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoVectorial2D U;
VF::TCampoEscalar2D p,
                    k,
                    ω,
                    νt;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<VF::TDirichlet>("inlet", U0);
k.DefCC<VF::TDirichlet>("inlet", k0);
ω.DefCC<VF::TDirichlet>("inlet", ω0);

U.DefCC<VF::TDirichlet>("airfoil");
k.DefCC<VF::TDirichlet>("airfoil");
ω.DefCC<TCCω>("airfoil");

U.DefCC<VF::TSimetria>("top");
U.DefCC<VF::TSimetria>("bottom");

p.DefCC<VF::TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;

solve(lap(p) == div(U), p);

U -= grad(p);

p = 0.0;
k = k0;
ω = ω0;

while (true)
  {
// ------------------------------------------------------------------------------- Campos calculados

  νt = k / ω;

// ---------------------------------------------------------------------------------- Momento lineal

  VF::TSistema const UEc = div(U * U) - div((ν + νt) * grad(U)) == (grad(νt) & gradT(U)) - grad(p);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  VF::TCampo const pOld = p;

  U = VF::TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

  for (std::size_t i = 0u; i < NCorNoOrto; ++i)
    solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 0.001 * sum(mag(p)))
    break;

  p = lerp(p, pOld, α);

// ----------------------------------------------------------------------------- Transporte de k y ω

  auto const G = νt * ((grad(U) + gradT(U)) && grad(U));

  solve(div(U * k) - div((ν + σk * νt) * grad(k)) + Cμ * ω * (+k) ==          G, k, α);

  solve(div(U * ω) - div((ν + σω * νt) * grad(ω)) + β  * ω * (+ω) == γ / νt * G, ω, α);
  }

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");
k.Write(ofs, "k");
ω.Write(ofs, "omega");
νt.Write(ofs, "nut");

return 0;
}
