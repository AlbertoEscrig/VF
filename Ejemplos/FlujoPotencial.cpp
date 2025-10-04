// =================================================================================================
// ============================================================================== FlujoPotencial.cpp
//
// Flujo estacionario incompresible invíscido irrotacional en torno a un cilindro
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

// -------------------------------------------------------------------------------------- Constantes

constexpr VF::TVector2D U0 = {1.0, 0.0};

constexpr std::size_t NCorNoOrto = 3u;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla2D::Read("Cylinder.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoVectorial2D U;
VF::TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<VF::TDirichlet>("left", U0);
U.DefCC<VF::TSimetria>("up");
U.DefCC<VF::TSimetria>("down");
U.DefCC<VF::TSimetria>("cylinder");

p.DefCC<VF::TDirichlet>("right");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;

// ---------------------------------------------------------------------- Corrección de la velocidad

do
  {
  for (std::size_t i = 0u; i < NCorNoOrto; ++i)
    solve(lap(p) == div(U), p);

  U -= grad(p);
  }
while (sum(mag(grad(p))) > 1e-3 * sum(mag(U)));

// ------------------------------------------------------------------------------ Campo de presiones

VF::TCampo const divUU = (U * U) / (U & U) & div(U * U);

p = 0.0;

for (std::size_t i = 0u; i < NCorNoOrto; ++i)
  solve(lap(p) == -div(divUU), p);

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");

return 0;
}
