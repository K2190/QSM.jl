#pragma once

#include <stddef.h>

extern void coarsen_grid(uint8_t* G2,
    const uint8_t* G, const size_t* sz2, const size_t* sz);
