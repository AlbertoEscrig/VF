// =================================================================================================
// ====================================================================================== SIMPLE.cpp
//
// SIMPLE (Semi Implicit Method for Pressure Linked Equations)
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector2D U0 = {1.0, 0.0};

constexpr double ν = 0.01;

constexpr double α = 0.7;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("Cavity.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("moving wall", U0);
U.DefCC<TDirichlet>("fixed walls");

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

  TSistema pEc = div((1.0 / UEc.aP) * grad(p)) == div(U);

  pEc.DefRef(0u, 0.0);

  solve(pEc, p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 0.1)
    break;

  p = lerp(p, pOld, α);
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");

return 0;
}
