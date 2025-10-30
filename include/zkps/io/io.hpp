#pragma once
#include "zkps/types.hpp"


namespace zkps::io {


void read_w(uint256* ws, data_stream* w_s1, data_stream* w_s2, data_stream* w_s3, data_stream* w_s4);
void read_s(uint256* scalars, s_stream* s1);
void read_coeffs(uint256* coeffs, data_stream* c_s, uint32 step);


void subcore(data_stream* s1, data_stream* s2, data_stream* s3);
void mulcore(data_stream* s1, data_stream* s2, data_stream* s3);


void writeback(uint256* nums, data_stream* n_s, uint32 step);


void lazy_read(uint256* p_1,uint256* p_2,uint256* p_3,uint256* p_4,uint256* p_5,uint256* p_6,
    data_stream0* x_s,data_stream0* y_s,data_stream0* z_s,data_stream0* t_s);


void lazy_write(uint256* p_1,uint256* p_2,uint256* p_3,uint256* p_4,uint256* p_5,uint256* p_6,
    data_stream0* x_s,data_stream0* y_s,data_stream0* z_s,data_stream0* t_s);
void readnum(uint256* num, data_stream* s, uint32 step);

// High‑level loader that wires all streams for NTT/MSM pipelines
void data_loader_dual(
    uint32 mode, uint32 step,
    uint256* mem0, uint256* mem1, uint256* mem2, uint256* mem3,
    uint256* mem4, uint256* mem5, uint256* mem6, uint256* mem7,
    uint256* mem8, uint256* mem9, uint256* mem10, uint256* mem11,
    uint256* mem12, uint256* mem13, uint256* mem14, uint256* mem15,
    uint256* mem16,
    data_stream* w_s1, data_stream* w_s2, data_stream* w_s3, data_stream* w_s4,
    data_stream* left0, data_stream* right0,
    data_stream* left1, data_stream* right1,
    data_stream* left2, data_stream* right2,
    data_stream* left3, data_stream* right3,
    data_stream* ca0, data_stream* cb0,
    data_stream* ca1, data_stream* cb1,
    data_stream* ca2, data_stream* cb2,
    data_stream* ca3, data_stream* cb3,
    s_stream* s1,
    data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s,
    data_stream0* xw_s, data_stream0* yw_s, data_stream0* zw_s, data_stream0* tw_s);


} // namespace zkps::io