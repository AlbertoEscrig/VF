// =================================================================================================
// =================================================================================== Catenoide.cpp
//
// Superficie de área mínima entre dos circunferencias concéntricas a distinta altura
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla2D::Read("Catenoide.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoEscalar2D z;

// ------------------------------------------------------------------------- Condiciones de contorno

z.DefCC<VF::TDirichlet>("inner", 0.0);
z.DefCC<VF::TDirichlet>("outer", 0.5);
z.DefCC<VF::TSimetria>("symmetry");

// ---------------------------------------------------------------------------------------- Solución

while (true)
  {
  VF::TCampo const zOld = z;

  VF::TCampo const γ = 1.0 / sqrt(1.0 + (grad(z) & grad(z)));

  solve(div(γ * grad(z)) == 0.0, z);

  if (sum(mag(z - zOld)) < 1e-5)
    break;
  }

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TMalla2D::Write(ofs);
z.Write(ofs, "z");

return 0;
}
