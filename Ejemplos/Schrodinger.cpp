// =================================================================================================
// ================================================================================= Schrodinger.cpp
//
// Orbital 3dz² del átomo de hidrógeno
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr int n = 3;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla3D::Read("Hidrogeno.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoEscalar3D Ψ;

// ------------------------------------------------------------------------- Condiciones de contorno

Ψ.DefCC<TDirichlet>("far field");
Ψ.DefCC<TSimetria>("symmetry");

// ------------------------------------------------------------------------------------- Expresiones

auto const dV = 8.0 * TMalla3D::V();   // Se simula un octante
auto const r  = mag(TMalla3D::C());    // Distancia al núcleo
auto const HΨ = -lap(Ψ) / 2.0 - Ψ / r; // Hamiltoniano

// --------------------------------------------------------------------------- Condiciones iniciales

Ψ = 1.0 / sqrt(sum(dV));

double E = sum(Ψ * HΨ * dV);

while (true)
  {
  double const EOld = E;

// ---------------------------------------------------------------------------------------- Solución

  solve(-lap(Ψ) / 2.0 - (1.0 / r - 1.0 / (2.0 * n * n)) * Sp(Ψ) == Ψ, Ψ);

  Ψ /= sqrt(sum(Ψ * Ψ * dV));

  E = sum(Ψ * HΨ * dV);

  if (abs(E - EOld) < 1e-7)
    break;
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla3D::Write(ofs);
Ψ.Write(ofs, "Psi");

return 0;
}
