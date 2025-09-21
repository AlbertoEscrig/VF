// =================================================================================================
// =============================================================================== TuboPeriodica.cpp
//
// Flujo laminar en un tubo con condiciones de contorno periódicas
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector3D U0 = {1.0, 0.0, 0.0};

constexpr double Δp = 10.0;

constexpr double ν = 0.01;

constexpr double α = 0.7;

constexpr TTensor<3u, 2u> T = {1.0,  0.0, 0.0,
                               0.0,  0.0, 1.0,
                               0.0, -1.0, 0.0};

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla3D::Read("TuboPeriodica.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial3D U;
TCampoEscalar3D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TPeriodica>("inlet");
U.DefCC<TPeriodica>("outlet");
U.DefCC<TDirichlet>("wall");
U.DefCC<TPeriodica>("xy plane", T);
U.DefCC<TPeriodica>("xz plane", T.T());

p.DefCC<TDirichlet>("inlet", Δp);
p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;

while (true)
  {
// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc = div(U * U) - ν * lap(U) == -grad(p);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  U = TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

  TCampo const pOld = p;

  solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 1e-3)
    break;

  p = lerp(p, pOld, α);
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla3D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");

return 0;
}
