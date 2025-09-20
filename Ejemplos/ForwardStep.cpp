// =================================================================================================
// ================================================================================= ForwardStep.cpp
//
// Flujo supersónico en torno a un obstáculo
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double μ  = 18.1e-6,  // Pa s
                 R  = 0.714286, // J / kg K
                 Cv = 1.78571,  // J / kg K
                 λ  = 32.3e-6;  // W / m K

constexpr TVector2D U0 = {3.0, 0.0}; // Mach 3
constexpr double p0 = 1.0,
                 T0 = 1.0;

constexpr double dt      = 0.002,
                 dtWrite = 0.5,
                 tFin    = 5.0;

constexpr auto I = TTensor<2u, 2u>::I();

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("ForwardStep.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;

TCampoEscalar2D p,
                T,
                ρ;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
p.DefCC<TDirichlet>("inlet", p0);
T.DefCC<TDirichlet>("inlet", T0);
ρ.DefCC<TDirichlet>("inlet", p0 / (R * T0));
U.DefCC<TSimetria>("top");
U.DefCC<TSimetria>("bottom");
U.DefCC<TDirichlet>("obstacle");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;
p = p0;
T = T0;

for (double t = dt; t < tFin + dt; t += dt)
  {
  TCampo const UOld = U;
  TCampo const pOld = p;
  TCampo const TOld = T;

  while (true)
    {
// ------------------------------------------------------------------------------- Campos calculados

    ρ = p / (R * T);

// ---------------------------------------------------------------------------------- Momento lineal

    TSistema const UEc =
      ρ * Sp(U) / dt - div(ρ * U) * Sp(U) + div(ρ * U * U) - μ * lap(U)
      ==
      ρ * UOld / dt + μ * grad(div(U)) / 3.0 - grad(p);

    solve(UEc, U);

// ------------------------------------------------------------------------------------- Continuidad

    TCampo const pNew = p;

    U = TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

    TSistema const pEc =
      Sp(p) / (R * T) / dt + div((U / (R * T)) * p) - div((ρ / UEc.aP) * grad(p))
      ==
      pOld / (R * TOld) / dt - div(ρ * U);

    solve(pEc, p);

    U -= grad(p) / UEc.aP;

    if (sum(mag(p - pNew)) < 1e-5)
      break;

// ----------------------------------------------------------------------------------------- Energía

    auto const G = μ * ((grad(U) + gradT(U) - 2.0 / 3.0 * div(U) * I) && grad(U));

    TSistema const TEc =
      Cv * (ρ * Sp(T) / dt - div(ρ * U) * Sp(T) + div(ρ * U * T)) - λ * lap(T)
      ==
      ρ * Cv * TOld / dt + G - p * div(U);

    solve(TEc, T);
    }

// -------------------------------------------------------------------------------------- Resultados

  if (fmod(t, dtWrite) < dt)
    {
    static int i = 0;

    ofstream ofs(format("resu{:02d}.vtk", i++));

    TMalla2D::Write(ofs);
    U.Write(ofs, "U");
    p.Write(ofs, "p");
    T.Write(ofs, "T");
    }
  }

return 0;
}
