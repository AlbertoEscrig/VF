// =================================================================================================
// ====================================================================================== Ciclon.cpp
//
// Flujo en un ciclón
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

#include "kEpsilon.h"

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector3D U0 = {10.0, 0.0, 0.0};

constexpr double k0 = 0.375,
                 ε0 = 6.41;

constexpr double ν = 1.5e-5;

constexpr double α = 0.7;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla3D::Read("Ciclon.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial3D U;
TCampoEscalar3D p,
                k,
                ε,
                νt,
                νEf;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
k.DefCC<TDirichlet>("inlet", k0);
ε.DefCC<TDirichlet>("inlet", ε0);

U.DefCC<TDirichlet>("wall");
k.DefCC<TkWallFunc>("wall");
ε.DefCC<TεWallFunc>("wall", k);
νEf.DefCC<TμEfWallFunc>("wall", k, ν);

p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

solve(lap(p) == div(U), p);

U -= grad(p);

p = 0.0;
k = k0;
ε = ε0;

while (true)
  {
// ---------------------------------------------------------------------------------- Momento lineal

  νt = Cμ * k * k / ε;

  νEf = ν + νt;

  TSistema const UEc = div(U * U) - div(νEf * grad(U)) == (grad(νt) & gradT(U)) - grad(p);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  U = TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

  TCampo const pOld = p;

  solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 0.001 * sum(mag(p)))
    break;

  p = lerp(p, pOld, α);

// ----------------------------------------------------------------------------- Transporte de k y ε

  TCampo const ω = ε / k;

  auto const G = νt * ((grad(U) + gradT(U)) && grad(U));

  solve(div(U * k) - div((ν + νt / σk) * grad(k)) +      ω * Sp(k) ==          G, k, α);

  solve(div(U * ε) - div((ν + νt / σε) * grad(ε)) + C2 * ω * Sp(ε) == C1 * ω * G, ε, α);
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla3D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");
k.Write(ofs, "k");
ε.Write(ofs, "epsilon");
νt.Write(ofs, "nut");

return 0;
}
