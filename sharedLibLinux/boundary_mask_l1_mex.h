#pragma once

#include <inttypes.h>
#include <stdint.h>
#include <stddef.h>

extern void boundary_mask_l1(uint8_t * B, const uint8_t * G,
	                                              const size_t * sz, const uint8_t l1);
