// =================================================================================================
// ==================================================================================== Hartmann.cpp
//
// Flujo magnetohidrodinámico
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector2D U0 = {1.0, 0.0};
constexpr TVector2D B0 = {0.0, 20.0};

constexpr double ν = 1.0; // m² s⁻¹
constexpr double ρ = 1.0; // kg m⁻³
constexpr double σ = 1.0; // Ω⁻¹ m⁻¹
constexpr double μ = 1.0; // H m⁻¹

constexpr double α = 0.7;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("Hartmann.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U,
                  B;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("inlet", U0);
U.DefCC<TDirichlet>("walls");
B.DefCC<TDirichlet>("walls", B0);
p.DefCC<TDirichlet>("outlet");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;
B = B0;

while (true)
  {
// ---------------------------------------------------------------------------------- Momento lineal

  TSistema const UEc =
    div(U * U) - ν * lap(U) == div(B * B) / (μ * ρ) + grad(B & B) / (2.0 * μ * ρ) - grad(p);

  solve(UEc, U, α);

// ------------------------------------------------------------------------------------- Continuidad

  U = TCampo((UEc.b + grad(p) - UEc.ΣaN(U)) / UEc.aP);

  TCampo const pOld = p;

  solve(div((1.0 / UEc.aP) * grad(p)) == div(U), p);

  U -= grad(p) / UEc.aP;

  if (sum(mag(p - pOld)) < 0.1)
    break;

  p = lerp(p, pOld, α);

// --------------------------------------------------------------------- Densidad de flujo magnético

  solve(div(U * B) - lap(B) / (μ * σ) == div(B * U), B);
  }

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla2D::Write(ofs);
U.Write(ofs, "U");
B.Write(ofs, "B");
p.Write(ofs, "p");

return 0;
}
