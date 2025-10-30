#pragma once
#include <ap_int.h>
#include <hls_stream.h>
#include "zkps/config.hpp"


namespace zkps {


using uint1 = ap_uint<1>;
using uint32 = ap_uint<32>;
using inddata = ap_uint<WINDOW_SIZE>;
using uint128 = ap_uint<128>;
using uint256 = ap_uint<256>;
using uint384 = ap_uint<384>;
using uint512 = ap_uint<512>;
using uint768 = ap_uint<768>;


struct epoint {
// Jacobian point representation
uint384 x = 0;
uint384 y = 1;
uint384 z = 0;
uint384 t = 0; // optional
};


using ep_stream = hls::stream<epoint>;
using data_stream = hls::stream<uint256>;
using data_stream0 = hls::stream<uint384>;
using s_stream = hls::stream<inddata>;


} // namespace zkps