// =================================================================================================
// ==================================================================================== OneEqLES.cpp
//
// Large Eddy Simulation (LES) -- Modelo de turbulencia de Yoshizawa y Horiuti
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double dt      = 1e-5,
                 dtWrite = 1e-2,
                 tFin    = 0.2;

constexpr TVector2D U0 = {10.0, 0.0};

constexpr double k0 = 2e-5;

constexpr double ν = 1.8e-5;

constexpr double Ck = 0.094,
                 Cε = 1.048;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("BackwardStep.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p,
                k;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
k.DefCC<TDirichlet>("inlet", k0);
U.DefCC<TDirichlet>("wall");
k.DefCC<TDirichlet>("wall");
p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;
k = k0;

// ------------------------------------------------------------------------------- Campos calculados

auto const Δ = sqrt(TMalla2D::V()); // 2D
//auto const Δ = cbrt(TMalla3D::V()); // 3D

for (double t = dt; t < tFin + dt; t += dt)
  {
  TCampo const νt = Ck * sqrt(k) * Δ;

  TCampo const gradp = grad(p);

// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc =
    d(U) / dt + div(U * U) - div((ν + νt) * grad(U)) == (grad(νt) & gradT(U)) - gradp;

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

// --------------------------------------------------------------------- Energía cinética turbulenta

  auto const G = νt * ((grad(U) + gradT(U)) && grad(U));

  solve(d(k) / dt + div(U * k) - div((ν + νt) * grad(k)) + (Cε * sqrt(k) / Δ) * Sp(k) == G, k);

  k = max(k, 0.0);

// -------------------------------------------------------------------------------------- Resultados

  if (fmod(t, dtWrite) < dt)
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
