#pragma once
#include "zkps/types.hpp"


namespace zkps::core {


// Stage 1 of NTT: butterfly with twiddle factors
void ntt1(
    data_stream* left0, data_stream* right0,
    data_stream* left1, data_stream* left2,
    data_stream* right1, data_stream* right2,
    data_stream* w_s, // twiddle stream
    data_stream* outa_s, data_stream* outb_s);


// Stage 2 of NTT: routing results back to left/right banks
void ntt2(
    data_stream* outa_s, data_stream* outb_s,
    data_stream* left1, data_stream* left2,
    data_stream* right1, data_stream* right2,
    data_stream* c1, data_stream* c2);


}