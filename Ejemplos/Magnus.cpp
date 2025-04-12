// =================================================================================================
// ====================================================================================== Magnus.cpp
//
// Efecto Magnus
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// --------------------------------------------------------------------------------------- FnRotWall

TVector<2u>
FnRotWall(TCara<2u> const &Cara)
{
constexpr double vr = 10.0 * 2.0 * numbers::pi * 10e-2; // 10 Hz
TVector<2u> const nf = Cara.nf();

return {vr * nf[1u], -vr * nf[0u]};
}

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector2D U0 = {0.0, 1.0};

constexpr double ν = 0.01;

constexpr double α = 0.7;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("Magnus.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
U.DefCC<TFuncion>("ball", &FnRotWall);
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
