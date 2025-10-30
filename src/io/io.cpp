// src/io/io.cpp
#include "zkps/io/io.hpp"
#include "zkps/constants.hpp"     // bls377_p / bls377_i / etc.
#include "zkps/math/field.hpp"   // SUB / MUL helpers if you prefer, but we call inline versions here
#include "zkps/config.hpp"

namespace zkps::io {

using namespace zkps;
using namespace zkps::consts;
using namespace zkps::math;


void read_w(uint256* ws, data_stream* w_s1, data_stream* w_s2, data_stream* w_s3, data_stream* w_s4) {
read_w_label1:
    for (int i = 0; i < R * LOGR / 2; i++) {
    #pragma HLS PIPELINE
        uint256 tempw = ws[i];
        w_s1->write(tempw);
        w_s2->write(tempw);
        w_s3->write(tempw);
        w_s4->write(tempw);
    }
}


void read_s(uint256* scalars, s_stream* s1) {
    const uint256 ONE = 1;
    for (uint32 i = 0; i < S; i++) {
    #pragma HLS PIPELINE II=57
        uint256 temps = scalars[i];
        inddata s_ary[WINDOW_NUM];
        for (int j = 0; j < WINDOW_NUM - 1; j++) {
        #pragma HLS PIPELINE off
            ap_uint<WINDOW_SIZE> temp = temps((j + 1) * WINDOW_SIZE - 1, j * WINDOW_SIZE);
            if (temp[WINDOW_SIZE - 1] == 1) {
                temps += ONE << ((j + 1) * WINDOW_SIZE);
            }
            s_ary[j] = temp(WINDOW_SIZE - 1, 0);
        }
        s_ary[WINDOW_NUM - 1] = temps(255, (WINDOW_NUM - 1) * WINDOW_SIZE);
        for (int k = 0; k < WINDOW_NUM; k++) {
            s1->write(s_ary[k]);
        }
    }
}


void read_coeffs(uint256* coeffs, data_stream* c_s, uint32 step) {
    uint32 startadd = (step * R) >> 1;
read_coeffs_label1:
    for (int i = 0; i < R / 2; i++) {
    #pragma HLS PIPELINE
        uint256 tempc = coeffs[startadd + i];
        c_s->write(tempc);
    }
}

void subcore(data_stream* s1, data_stream* s2, data_stream* s3) {
subcore_label0:
    for (int i = 0; i < R; i++) {
    #pragma HLS PIPELINE
        uint256 t1 = s1->read();
        uint256 t2 = s2->read();
        uint256 t3 = SUB(t1, t2, bls377_p);
        s3->write(t3);
    }
}

void mulcore(data_stream* s1, data_stream* s2, data_stream* s3) {
mulcore_label0:
    for (int i = 0; i < R; i++) {
    #pragma HLS PIPELINE
        uint256 t1 = s1->read();
        uint256 t2 = s2->read();
        uint256 t3 = MUL(t1, t2, bls377_p, bls377_i);
        s3->write(t3);
    }
}

void writeback(uint256* nums, data_stream* n_s, uint32 step) {
    uint32 startadd = (step * R) >> 1;
    for (int i = 0; i < R / 2; i++) {
    #pragma HLS PIPELINE
        uint256 tep = n_s->read();
        nums[startadd + i] = tep;
    }
}

void lazy_read(uint256* p_1, uint256* p_2, uint256* p_3, uint256* p_4, uint256* p_5, uint256* p_6,
               data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s) {
    for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE
        uint384 xt = (uint384(p_1[i]) << 128) + p_5[i](127, 0);
        uint384 yt = (uint384(p_2[i]) << 128) + p_5[i](255, 128);
        uint384 zt = (uint384(p_3[i]) << 128) + p_6[i](127, 0);
        uint384 tt = (uint384(p_4[i]) << 128) + p_6[i](255, 128);
        x_s->write(xt);
        y_s->write(yt);
        z_s->write(zt);
        t_s->write(tt);
    }
}

void lazy_write(uint256* p_1, uint256* p_2, uint256* p_3, uint256* p_4, uint256* p_5, uint256* p_6,
                data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s) {
    for (int i = 0; i < 32; i++) {
    #pragma HLS PIPELINE
        uint384 xt = x_s->read();
        uint384 yt = y_s->read();
        uint384 zt = z_s->read();
        uint384 tt = t_s->read();

        uint256 p1 = xt(383, 128);
        uint256 p2 = yt(383, 128);
        uint256 p3 = zt(383, 128);
        uint256 p4 = tt(383, 128);
        uint256 p5 = uint256(xt(127, 0)) + (uint256(yt(127, 0)) << 128);
        uint256 p6 = uint256(zt(127, 0)) + (uint256(tt(127, 0)) << 128);

        p_1[i] = p1;  p_2[i] = p2;  p_3[i] = p3;
        p_4[i] = p4;  p_5[i] = p5;  p_6[i] = p6;
    }
}


void data_loader_dual(
    uint32 mode, uint32 step,
    uint256* mem0,  uint256* mem1,  uint256* mem2,  uint256* mem3,
    uint256* mem4,  uint256* mem5,  uint256* mem6,  uint256* mem7,
    uint256* mem8,  uint256* mem9,  uint256* mem10, uint256* mem11,
    uint256* mem12, uint256* mem13, uint256* mem14, uint256* mem15,
    uint256* mem16,
    data_stream* w_s1, data_stream* w_s2, data_stream* w_s3, data_stream* w_s4,
    data_stream* left0,  data_stream* right0,
    data_stream* left1,  data_stream* right1,
    data_stream* left2,  data_stream* right2,
    data_stream* left3,  data_stream* right3,
    data_stream* ca0, data_stream* cb0,
    data_stream* ca1, data_stream* cb1,
    data_stream* ca2, data_stream* cb2,
    data_stream* ca3, data_stream* cb3,
    s_stream*    s1,
    data_stream0* x_s,  data_stream0* y_s,  data_stream0* z_s,  data_stream0* t_s,
    data_stream0* xw_s, data_stream0* yw_s, data_stream0* zw_s, data_stream0* tw_s
) {
#pragma HLS DATAFLOW
    if (mode == 0) {
        // NTT branch
        read_w(mem0, w_s1, w_s2, w_s3, w_s4);

        read_coeffs(mem1, left0, step);
        read_coeffs(mem2, right0, step);

        subcore(left0, right0, left1);
        mulcore(left1, w_s1, right1);

        writeback(mem9,  ca0, step);
        writeback(mem10, cb0, step);
        writeback(mem11, ca1, step);
        writeback(mem12, cb1, step);
        writeback(mem13, ca2, step);
        writeback(mem14, cb2, step);
        writeback(mem15, ca3, step);
        writeback(mem16, cb3, step);
    } else {
        read_s(mem0, s1);
        lazy_read(mem1, mem2, mem3, mem4, mem5, mem6, x_s,  y_s,  z_s,  t_s);
        lazy_write(mem9, mem10, mem11, mem12, mem13, mem14, xw_s, yw_s, zw_s, tw_s);
    }
}
void readnum(uint256* num, data_stream* s, uint32 step) {
    // Stream R/2 numbers starting at (step*R)/2
    uint32 start = (step * R) >> 1;
readnum_loop:
    for (int i = 0; i < R/2; ++i) {
    #pragma HLS PIPELINE
        s->write(num[start + i]);
    }
}
} // namespace zkps::io
