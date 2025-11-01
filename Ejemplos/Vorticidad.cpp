// =================================================================================================
// ================================================================================== Vorticidad.cpp
//
// Interacción entre vórtices perpendiculares
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

// -------------------------------------------------------------------------------------- Constantes

constexpr double ν  = 0.001,
                 Γ  = 1.0,
                 r0 = 0.1,
                 d0 = 0.3;

constexpr VF::TVector3D i = {1.0, 0.0, 0.0},
                        j = {0.0, 1.0, 0.0},
                        k = {0.0, 0.0, 1.0};

constexpr double Δt      = 0.005,
                 ΔtWrite = 0.1,
                 tFin    = 4.0;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla3D::Read("CuboPeriodico.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoVectorial3D U,
                      ω;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<VF::TPeriodica>("xy plane");
U.DefCC<VF::TPeriodica>("xz plane");
U.DefCC<VF::TPeriodica>("yz plane");

ω.DefCC<VF::TPeriodica>("xy plane");
ω.DefCC<VF::TPeriodica>("xz plane");
ω.DefCC<VF::TPeriodica>("yz plane");

// --------------------------------------------------------------------------- Condiciones iniciales

auto const C = VF::TMalla3D::C();

auto const ω1 = exp(-(pow<2u>(C & j) + pow<2u>((C & k) - d0)) / (r0 * r0)) * i;
auto const ω2 = exp(-(pow<2u>(C & i) + pow<2u>((C & k) + d0)) / (r0 * r0)) * j;

ω = Γ / (std::numbers::pi * r0 * r0) * (ω1 + ω2);

solve(lap(U) == -rot(ω), U);

for (double t = Δt; t < tFin + Δt; t += Δt)
  {
// ---------------------------------------------------------------------------------------- Solución

  solve(Δ(ω) / Δt + div(U * ω) - ν * lap(ω) == (ω & grad(U)), ω);

  solve(lap(U) == -rot(ω), U);

// -------------------------------------------------------------------------------------- Resultados

  if (std::fmod(t, ΔtWrite) < Δt)
    {
    static int i = 0;

    std::ofstream ofs(std::format("resu{:02d}.vtk", i++));

    VF::TMalla3D::Write(ofs);
    U.Write(ofs, "U");
    ω.Write(ofs, "omega");
    }
  }

return 0;
}
