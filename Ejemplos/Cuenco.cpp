// =================================================================================================
// ====================================================================================== Cuenco.cpp
//
// Transmisión de calor en un cuenco parabólico (z = a * (x² + y²))
//
// ================================================================= Copyright © 2025 Alberto Escrig
// =================================================================================================

import VF;
import std;


// -------------------------------------------------------------------------------------- Constantes

constexpr double a  = 10.0,
                 κ2 = 4.0 * a * a;

constexpr double Thot  = 10.0,
                 Tcold = 0.0;

constexpr auto I = VF::TTensor<2u, 2u>::I();

constexpr std::size_t NCorNoOrto = 5u;

int main()
{
// ------------------------------------------------------------------------------------------- Malla

VF::TMalla2D::Read("Cuenco.msh");

// ------------------------------------------------------------------------------------------ Campos

VF::TCampoEscalar2D T;

// ------------------------------------------------------------------------- Condiciones de contorno

T.DefCC<VF::TDirichlet>("hot", Thot);
T.DefCC<VF::TDirichlet>("cold", Tcold);
T.DefCC<VF::TSimetria>("symmetry");

// ---------------------------------------------------------------------------------------- Solución

auto const r = VF::TMalla2D::C();

auto const K = ((1.0 + κ2 * (r & r)) * I - κ2 * (r * r)) / sqrt(1.0 + κ2 * (r & r));

for (std::size_t i = 0u; i < NCorNoOrto; ++i)
  solve(div(K & grad(T)) == 0.0, T);

// -------------------------------------------------------------------------------------- Resultados

std::ofstream ofs("resu.vtk");

VF::TCampo const q = -K & grad(T);

VF::TMalla2D::Write(ofs);
T.Write(ofs, "T");
q.Write(ofs, "q");

return 0;
}
