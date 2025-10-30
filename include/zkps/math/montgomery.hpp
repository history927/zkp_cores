#pragma once
#include "zkps/types.hpp"


namespace zkps::math {


uint256 montgomery_reduce_256(const uint512& t, const uint256& m, const uint256& inv);
uint384 montgomery_reduce_384(const uint768& t, const uint384& m, const uint384& inv);


inline uint256 montgomery_reduce(const uint512& t, const uint256 m, const uint256 inv) {
    #pragma HLS INLINE
    return montgomery_reduce_256(t, m, inv);
}
inline uint384 montgomery_reduce0(const uint768& t, const uint384& m, const uint384& inv) {
    #pragma HLS INLINE
    return montgomery_reduce_384(t, m, inv);
}


}