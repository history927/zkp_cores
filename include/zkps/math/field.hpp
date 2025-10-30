#pragma once
#include "zkps/types.hpp"
#include "zkps/math/montgomery.hpp"


namespace zkps::math {


// Modular add/sub for 256 and 384 bits
inline uint256 add_mod(uint256 a, uint256 b, const uint256& q) {
    #pragma HLS INLINE
    uint256 t = a + b; return (t >= q) ? (t - q) : t;
}
inline uint256 sub_mod(uint256 a, uint256 b, const uint256& q) {
    #pragma HLS INLINE
    return (a >= b) ? (a - b) : (q - (b - a));
}
inline uint384 add_mod(uint384 a, uint384 b, const uint384& q) {
    #pragma HLS INLINE
    uint384 t = a + b; return (t >= q) ? (t - q) : t;
}
inline uint384 sub_mod(uint384 a, uint384 b, const uint384& q) {
    #pragma HLS INLINE
    return (a >= b) ? (a - b) : (q - (b - a));
}


// Montgomery mul for 256 and 384 bits
inline uint256 mul_mont(uint256 a, uint256 b, const uint256& m, const uint256& inv) {
    #pragma HLS INLINE
    return montgomery_reduce_256(karatsuba_256(a, b), m, inv);
}
inline uint384 mul_mont(uint384 a, uint384 b, const uint384& m, const uint384& inv) {
    #pragma HLS INLINE
    return montgomery_reduce_384(karatsuba_384(a, b), m, inv);
}


// === Compatibility layer (preserve original function names) ===
inline uint256 ADD(uint256 a, uint256 b, uint256 q) { return add_mod(a, b, q); }
inline uint256 SUB(uint256 a, uint256 b, uint256 q) { return sub_mod(a, b, q); }
inline uint384 ADD0(uint384 a, uint384 b, const uint384 q) { return add_mod(a, b, q); }
inline uint384 SUB0(uint384 a, uint384 b, const uint384 q) { return sub_mod(a, b, q); }
inline uint256 MUL(uint256 a, uint256 b, const uint256 m, const uint256 inv) { return mul_mont(a, b, m, inv); }
inline uint384 MUL0(uint384 a, uint384 b, const uint384 q, const uint384 inv) { return mul_mont(a, b, q, inv); }


} // namespace zkps::math