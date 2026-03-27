/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef __MACHINE_VECTORS__H
#define __MACHINE_VECTORS__H

#include <flint/flint.h>

#include "pml.h"

/* Disable all machine vectors code to avoid compilation issues */
#undef PML_HAVE_MACHINE_VECTORS
#undef PML_HAVE_AVX2
#undef PML_HAVE_AVX512

#endif /* ifndef __MACHINE_VECTORS__H */
