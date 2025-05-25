# VF

This is a tiny C++ module library for finite-volume calculations
based on expression templates and lazy evaluation.

There are usage examples in the [Ejemplos](Ejemplos/) folder, which
should be understandable by those familiar with the finite-volume
framework applied to transport phenomena problems.

The library tries to mimic the mathematical language. For example,
solving the momentum equation:
```math
\nabla·\mathbf{U}\mathbf{U}-\nu\nabla·\nabla\mathbf{U}=-\nabla p
```
can be written:
```c++
solve(div(U * U) - ν * lap(U) == -grad(p), U);
```
The formalism is that every term in the left-hand side is evaluated
implicitly, whereas the right-hand-side terms are evaluated explicitly.

The library requires the [Gmsh](https://gmsh.info) C++ API to read the
input mesh. On Debian-based systems, it is part of the `libgmsh-dev`
package.

This library is distributed under the terms of the GNU General Public
License (GPL). See [LICENSE.txt](LICENSE.txt).

Copyright © 2025 Alberto Escrig

> [!IMPORTANT]
> This library is free software: you can redistribute it and/or modify
> it under the terms of the GNU General Public License as published by
> the Free Software Foundation, either version 3 of the License, or
> (at your option) any later version.
>
> This library is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU General Public License for more details.
