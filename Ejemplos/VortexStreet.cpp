// =================================================================================================
// ================================================================================ VortexStreet.cpp
//
// Flujo en torno a un cilindro en régimen no estacionario
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double dt      = 2e-5,
                 dtWrite = 1e-2,
                 tFin    = 0.6;

constexpr VF::TVector2D U0 = {1.0, 0.0},
                        δU = {0.0, 0.01};

constexpr double ν = 1e-5;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla2D::Read("VortexStreet.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoVectorial2D U;
VF::TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<VF::TDirichlet>("inlet", U0);
U.DefCC<VF::TDirichlet>("cylinder");
U.DefCC<VF::TSimetria>("up");
U.DefCC<VF::TSimetria>("down");
p.DefCC<VF::TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0 + δU;

for (double t = dt; t < tFin + dt; t += dt)
  {
// ------------------------------------------------------------------------------- Campos calculados

  VF::TCampo const gradp = grad(p);

// ---------------------------------------------------------------------------------- Momento lineal

  VF::TSistema const UEc = d(U) / dt + div(U * U) - ν * lap(U) == -gradp;

  solve(UEc, U);

// ------------------------------------------------------------------------------------- Continuidad

  while (true)
    {
    VF::TCampo const pOld = p;

    U = VF::TCampo((UEc.b + gradp - UEc.ΣaN(U)) / UEc.aP);

    solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

    U -= grad(p) / UEc.aP;

    if (sum(mag(p - pOld)) < 0.001 * sum(mag(p)))
      break;
    }

// -------------------------------------------------------------------------------------- Resultados

  if (std::fmod(t, dtWrite) < dt)
    {
    static int i = 0;

    std::ofstream ofs(std::format("resu{:02d}.vtk", i++));

    VF::TMalla2D::Write(ofs);
    U.Write(ofs, "U");
    p.Write(ofs, "p");
    }
  }

return 0;
}
