// =================================================================================================
// ================================================================================= Smagorinsky.cpp
//
// Modelo de turbulencia de Smagorinsky y Lilly
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double Δt      = 1e-5,
                 ΔtWrite = 1e-2,
                 tFin    = 0.2;

constexpr TVector2D U0 = {10.0, 0.0};

constexpr double ν = 1.5e-5;

constexpr double Cs = 0.15;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("BackwardStep.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
U.DefCC<TDirichlet>("wall");
p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;

for (double t = Δt; t < tFin + Δt; t += Δt)
  {
// ------------------------------------------------------------------------------- Campos calculados

  TCampo const νt = Cs * Cs * TMalla2D::V() * sqrt(2.0) * mag(grad(U) + gradT(U)); // 2D
  //TCampo const νt = pow<2u>(Cs * cbrt(TMalla3D::V())) * sqrt(2.0) * mag(grad(U) + gradT(U)); // 3D

  TCampo const gradp = grad(p);

// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc =
    Δ(U) / Δt + div(U * U) - div((ν + νt) * grad(U)) == (grad(νt) & gradT(U)) - gradp;

  solve(UEc, U);

// ------------------------------------------------------------------------------------- Continuidad

  while (true)
    {
    TCampo const pOld = p;

    U = TCampo((UEc.b + gradp - UEc.ΣaN(U)) / UEc.aP);

    solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

    U -= grad(p) / UEc.aP;

    if (sum(mag(p - pOld)) < 0.1)
      break;
    }

// -------------------------------------------------------------------------------------- Resultados

  if (fmod(t, ΔtWrite) < Δt)
    {
    static int i = 0;

    ofstream ofs(format("resu{:02d}.vtk", i++));

    TMalla2D::Write(ofs);
    U.Write(ofs, "U");
    p.Write(ofs, "p");
    νt.Write(ofs, "nut");
    }
  }

return 0;
}
