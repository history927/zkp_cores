#pragma once
#include "zkps/types.hpp"


namespace zkps::core {

void unipadd(
    data_stream0* x1_s, data_stream0* y1_s, data_stream0* z1_s, data_stream0* t1_s,
    data_stream0* x2_s, data_stream0* y2_s, data_stream0* z2_s, data_stream0* t2_s,
    data_stream0* xr_s, data_stream0* yr_s, data_stream0* zr_s, data_stream0* tr_s);
} // namespace zkps::core

// Three-input modular multiply pipeline (384-bit domain)
void threemul(data_stream0* in1, data_stream0* in2, data_stream0* out);


// Point-add pipeline stages (rename or extend if you use different models)
void padds1(data_stream0* x1_s, data_stream0* y1_s, data_stream0* t1_s,
    data_stream0* x2_s, data_stream0* y2_s, data_stream0* t2_s,
    data_stream0* out1, data_stream0* out2);


void padds2(data_stream0* m1_s, data_stream0* z1_s, data_stream0* z2_s,
    data_stream0* e_s, data_stream0* h_s,
    data_stream0* out1, data_stream0* out2);


void padds3(data_stream0* m2_s, data_stream0* e_s, data_stream0* h_s,
    data_stream0* out1, data_stream0* out2, data_stream0* t_s);


void padds4(data_stream0* m3_s, data_stream0* x_s, data_stream0* y_s, data_stream0* z_s);


// Bucketization helpers for MSM
void selector1(s_stream* s1, data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s,
    data_stream0* xr_s, data_stream0* yr_s, data_stream0* zr_s, data_stream0* tr_s);


void selector2(s_stream* s1, data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s,
    data_stream0* xr2_s, data_stream0* yr2_s, data_stream0* zr2_s, data_stream0* tr2_s);


void bucket_part(data_stream0* xr2_s, data_stream0* yr2_s, data_stream0* zr2_s, data_stream0* tr2_s,
    data_stream0* xw_s, data_stream0* yw_s, data_stream0* zw_s, data_stream0* tw_s);


} // namespace zkps::core