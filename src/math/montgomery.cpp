#include "zkps/math/montgomery.hpp"
#include "zkps/math/karatsuba.hpp"


namespace zkps::math {


uint256 montgomery_reduce_256(const uint512& t, const uint256& m, const uint256& inv) {
    #pragma HLS INLINE
    uint256 t_low = t(255, 0);
    uint512 q = karatsuba_256(t_low, inv);
    uint256 q0 = q(255, 0);
    uint512 x = t + karatsuba_256(q0, m); // may overflow; upper half is the candidate
    uint256 res = x(511, 256);
    if (res >= m) res = res - m;
    return res;
}


uint384 montgomery_reduce_384(const uint768& t, const uint384& m, const uint384& inv) {
    #pragma HLS INLINE
    uint384 t_low = t(383, 0);
    uint768 q = karatsuba_384(t_low, inv);
    uint384 q0 = q(383, 0);
    uint768 x = t + karatsuba_384(q0, m);
    uint384 res = x(767, 384);
    if (res >= m) res = res - m;
    return res;
}


} // namespace zkps::math