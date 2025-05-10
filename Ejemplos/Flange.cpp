// =================================================================================================
// ====================================================================================== Flange.cpp
//
// Conducción de calor en régimen no estacionario
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double dt      = 0.005,
                 dtWrite = 0.1,
                 tFin    = 3.0;

constexpr double α = 4e-5;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla3D::Read("Flange.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoEscalar3D T;

// ------------------------------------------------------------------------- Condiciones de contorno

T.DefCC<TDirichlet>("hot",  573.0);
T.DefCC<TDirichlet>("cold", 273.0);

// --------------------------------------------------------------------------- Condiciones iniciales

T = 273.0;

for (double t = dt; t < tFin + dt; t += dt)
  {
// ---------------------------------------------------------------------------------------- Solución

  solve(d(T) / dt - α * lap(T) == 0, T); // Euler implícito
  // solve(d(T) / dt - 0.5 * α * lap(T) == 0.5 * α * lap(T), T); // Crank-Nicolson
  // solve(d(T) / dt == α * lap(T), T); // Euler explícito

// -------------------------------------------------------------------------------------- Resultados

  if (fmod(t, dtWrite) < dt)
    {
    static int i = 0;

    ofstream ofs(format("resu{:02d}.vtk", i++));

    TMalla3D::Write(ofs);
    T.Write(ofs, "T");
    }
  }

return 0;
}
