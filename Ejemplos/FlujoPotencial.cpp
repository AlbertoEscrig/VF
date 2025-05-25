// =================================================================================================
// ============================================================================== FlujoPotencial.cpp
//
// Flujo estacionario incompresible invíscido irrotacional en torno a un cilindro
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

using namespace VF;
using namespace std;

// -------------------------------------------------------------------------------------- Constantes

constexpr TVector2D U0 = {1.0, 0.0};

constexpr size_t NCorNoOrto = 5u;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

TMalla2D::Read("Cylinder.msh");

// ------------------------------------------------------------------------------------------ Campos

TCampoVectorial2D U;
TCampoEscalar2D p;

// ------------------------------------------------------------------------- Condiciones de contorno

U.DefCC<TDirichlet>("left", U0);
U.DefCC<TSimetria>("up");
U.DefCC<TSimetria>("down");
U.DefCC<TSimetria>("cylinder");

p.DefCC<TDirichlet>("right");

// --------------------------------------------------------------------------- Condiciones iniciales

U = U0;

// ---------------------------------------------------------------------- Corrección de la velocidad

do
  {
  for (size_t i = 0u; i < NCorNoOrto; ++i)
    solve(lap(p) == div(U), p);

  U -= grad(p);
  }
while (sum(mag(p)) > 1e-2);

// ------------------------------------------------------------------------------ Campo de presiones

TCampo const divUU = (U * U) / (U & U) & div(U * U);

p = 0.0;

for (size_t i = 0u; i < NCorNoOrto; ++i)
  solve(lap(p) == -div(divUU), p);

// -------------------------------------------------------------------------------------- Resultados

ofstream ofs("resu.vtk");

TMalla2D::Write(ofs);
U.Write(ofs, "U");
p.Write(ofs, "p");

return 0;
}
