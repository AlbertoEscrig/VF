// =================================================================================================
// ===================================================================================== Airfoil.cpp
//
// Flujo en torno a un airfoil
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

#include "kEpsilon.h"

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double ν = 14e-6;

constexpr TVector2D U0 = {10.0, 0.0};

constexpr double k0 = 2.5e-8,
                 ε0 = 3.75e-8;

constexpr double α = 0.7;

constexpr size_t NCorNoOrto = 10u;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("naca5012.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p,
                k,
                ε,
                νt,
                νEf;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
k.DefCC<TDirichlet>("inlet", k0);
ε.DefCC<TDirichlet>("inlet", ε0);

U.DefCC<TSimetria>("top");
U.DefCC<TSimetria>("bottom");

U.DefCC<TDirichlet>("aerofoil");
k.DefCC<TkWallFunc>("aerofoil");
ε.DefCC<TεWallFunc>("aerofoil", k);
νEf.DefCC<TμEfWallFunc>("aerofoil", k, ν);

p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;
k = k0;
ε = ε0;

while (true)
  {
// ------------------------------------------------------------------------------- Campos calculados

  TCampo const ω = ε / k;

  νt = Cμ * k * k / ε;

  νEf = ν + νt;

// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc =
    div(U * U) - div(U) * Sp(U) - div(νEf * grad(U)) == (grad(νt) & gradT(U)) - grad(p);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  TCampo const pOld = p;

  U = TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

  for (size_t i = 0u; i < NCorNoOrto; ++i)
    solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 0.1)
    break;

  p = lerp(p, pOld, α);

// ----------------------------------------------------------------------------- Transporte de k y ε

  auto const G = νt * ((grad(U) + gradT(U)) && grad(U));

  solve(div(U * k) - div(U) * Sp(k) - div(νt * grad(k)) / σk +      ω * Sp(k) ==          G, k, α);

  solve(div(U * ε) - div(U) * Sp(ε) - div(νt * grad(ε)) / σε + C2 * ω * Sp(ε) == C1 * ω * G, ε, α);
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");
k.Write(ofs, "k");
ε.Write(ofs, "epsilon");
νt.Write(ofs, "nut");

return 0;
}
