// =================================================================================================
// ============================================================================== ForwardStepRho.cpp
//
// Flujo supersónico en torno a un obstáculo. Esquema basado en la densidad
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

constexpr double Δt      = 0.002,
                 ΔtWrite = 0.5,
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
T = T0;
ρ = p0 / (R * T0);

for (double t = Δt; t < tFin + Δt; t += Δt)
  {
  TCampo const UOld = U;
  TCampo const TOld = T;
  TCampo const ρOld = ρ;

  while (true)
    {
// ------------------------------------------------------------------------------- Campos calculados

    p = ρ * R * T;

// ---------------------------------------------------------------------------------- Momento lineal

    TSistema const UEc =
      ρ * Sp(U) / Δt - div(ρ * U) * Sp(U) + div(ρ * U * U) - μ * lap(U)
      ==
      ρ * UOld / Δt + μ * grad(div(U)) / 3.0 - grad(p);

    solve(UEc, U);

// ------------------------------------------------------------------------------------- Continuidad

    TSistema const ρEc = Sp(ρ) / Δt + div(U * ρ) == ρOld / Δt;

    solve(ρEc, ρ);

    if (sum(mag(p - ρ * R * T)) < 1e-3)
      break;

// ----------------------------------------------------------------------------------------- Energía

    auto const G = μ * ((grad(U) + gradT(U) - 2.0 / 3.0 * div(U) * I) && grad(U));

    TSistema const TEc =
      Cv * (ρ * Sp(T) / Δt - div(ρ * U) * Sp(T) + div(ρ * U * T)) - λ * lap(T)
      ==
      ρ * Cv * TOld / Δt + G - p * div(U);

    solve(TEc, T);
    }

// -------------------------------------------------------------------------------------- Resultados

  if (fmod(t, ΔtWrite) < Δt)
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
