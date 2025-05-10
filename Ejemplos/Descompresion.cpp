// =================================================================================================
// =============================================================================== Descompresion.cpp
//
// Descompresión de un tanque de agua presurizado
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector2D U0 = {};

constexpr double p0 = 100e5;     // Pa

constexpr double μ    = 1e-3,    // Pa s
                 pRef = 1e5,     // Pa
                 ρRef = 1e3,     // kg / m³
                 ψ    = 4.54e-7; // s² / m²

constexpr double dt      = 0.1e-6,
                 dtWrite =  10e-6,
                 tFin    = 250e-6;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("DecompressionTank.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

p.DefCC<TDirichlet>("nozzle");
U.DefCC<TDirichlet>("wall");
U.DefCC<TSimetria>("axis");

// --------------------------------------------------------------------------- Condiciones iniciales

p = p0;
U = U0;

for (double t = dt; t < tFin + dt; t += dt)
  {
  TCampo const UOld = U;
  TCampo const pOld = p;

  while (true)
    {
// ------------------------------------------------------------------------------- Campos calculados

    TCampo const ρ = ρRef + ψ * (p - pRef);

// ---------------------------------------------------------------------------------- Momento lineal

    TSistema const UEc =
      ρ * Sp(U) / dt - div(ρ * U) * Sp(U) + div(ρ * U * U) - μ * lap(U) == ρ * UOld / dt - grad(p);

    solve(UEc, U);

// ------------------------------------------------------------------------------------- Continuidad

    TCampo const pNew = p;

    U = TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

    TSistema const pEc =
      ψ * Sp(p) / dt + ψ * div(U * p) - div((ρ / UEc.aP) * grad(p)) == ψ * pOld / dt - div(ρ * U);

    solve(pEc, p);

    U -= grad(p) / UEc.aP;

    if (sum(mag(p - pNew)) < 0.1)
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
