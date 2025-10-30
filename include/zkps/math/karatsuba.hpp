#pragma once
#include "zkps/types.hpp"


namespace zkps::math {


// 128×128 → 256 (basic Karatsuba)
uint256 karatsuba_128(uint128 x, uint128 y);


// 256×256 → 512 (via 2-way split to 128)
uint512 karatsuba_256(uint256 x, uint256 y);


// 384×384 → 768 (via 3×128 split)
uint768 karatsuba_384(uint384 x, uint384 y);


} // namespace zkps::math