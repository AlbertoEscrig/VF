// =================================================================================================
// ================================================================================ VortexStreet.cpp
//
// Flujo en torno a un cilindro en régimen no estacionario
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double dt      = 2.5e-5,
                 dtWrite = 1e-2,
                 tFin    = 0.6;

constexpr TVector2D U0 = {1.0, 0.0};

constexpr double ν = 1e-5;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("VortexStreet.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
U.DefCC<TDirichlet>("cylinder");
U.DefCC<TSimetria>("up");
U.DefCC<TSimetria>("down");
p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;

for (double t = dt; t < tFin + dt; t += dt)
  {
// ------------------------------------------------------------------------------- Campos calculados

  TCampo const gradp = grad(p);

// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc = d(U) / dt + div(U * U) - ν * lap(U) == -gradp;

  solve(UEc, U);

// ------------------------------------------------------------------------------------- Continuidad

  while (true)
    {
    TCampo const pOld = p;

    U = TCampo((UEc.b + gradp - UEc.ΣaN(U)) / UEc.aP);

    solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

    U -= grad(p) / UEc.aP;

    if (sum(mag(p - pOld)) < 1e-3)
      break;
    }

// -------------------------------------------------------------------------------------- Resultados

  if (fmod(t, dtWrite) < dt)
    {
    static int i = 0;

    ofstream ofs(format("resu{:02d}.vtk", i++));

    TMalla2D::Write(ofs);
    U.Write(ofs, "U");
    p.Write(ofs, "p");
    }
  }

return 0;
}
