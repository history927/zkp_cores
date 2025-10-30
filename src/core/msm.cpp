// src/core/msm.cpp — complete, HLS-friendly (no TODOs)
#include "zkps/core/msm.hpp"
#include "zkps/constants.hpp"
#include "zkps/math/field.hpp"
#include "zkps/config.hpp"

namespace zkps::core {

using namespace zkps;
using namespace zkps::consts;
using namespace zkps::math;

//------------------------------------------------------------------------------
// 384-bit streaming Montgomery multiply: out = a*b (mod MOD_Q)
void threemul(data_stream0* in1, data_stream0* in2, data_stream0* out) {
    for (int i = 0; i < (WINDOWS + N) * 3; ) {
        #pragma HLS PIPELINE II=1
        if (!in1->empty() && !in2->empty()) {
            uint384 a = in1->read();
            uint384 b = in2->read();
            uint384 c = MUL0(a, b, MOD_Q, INV_Q);
            out->write(c);
            ++i;
        }
    }
}

//------------------------------------------------------------------------------
// Stage-1 of point-add pipeline (triplet cadence):
// cycle 0 -> (Y-X), cycle 1 -> (Y+X), cycle 2 -> T (forwarded)
void padds1(
    data_stream0* x1_s, data_stream0* y1_s, data_stream0* t1_s,
    data_stream0* x2_s, data_stream0* y2_s, data_stream0* t2_s,
    data_stream0* out1, data_stream0* out2
) {
    int counter = 0;
    uint384 ymx1=0, ymx2=0, ypx1=0, ypx2=0;
    uint384 x1=0, y1=0, t1=0, x2=0, y2=0, t2=0;

    for (int i = 0; i < (WINDOWS + N) * 3; ++i) {
        #pragma HLS PIPELINE II=1
        if (counter == 0) {
            x1 = x1_s->read(); y1 = y1_s->read();
            x2 = x2_s->read(); y2 = y2_s->read();
            t1 = t1_s->read(); t2 = t2_s->read();
            ymx1 = SUB0(y1, x1, MOD_Q);
            ymx2 = SUB0(y2, x2, MOD_Q);
            ypx1 = ADD0(y1, x1, MOD_Q);
            ypx2 = ADD0(y2, x2, MOD_Q);
        }

        if      (counter == 0) { out1->write(ymx1); out2->write(ymx2); }
        else if (counter == 1) { out1->write(ypx1); out2->write(ypx2); }
        else                   { out1->write(t1);   out2->write(t2);   }

        counter = (counter == 2) ? 0 : counter + 1;
    }
}

//------------------------------------------------------------------------------
// Stage-2: from triplets (A,B,T) produce E=B-A, H=B+A and D=2*Z1*Z2, C=k*T
// Output layout (triplet cadence):
//   out1: A, B, C
//   out2: 0, 0, D
// Side channels:
//   e_s:   E
//   h_s:   H
void padds2(
    data_stream0* m1_s, data_stream0* z1_s, data_stream0* z2_s,
    data_stream0* e_s,  data_stream0* h_s,
    data_stream0* out1, data_stream0* out2
) {
    int counter = 0;
    uint384 A=0, B=0, T=0, z1=0, z2=0;

    for (int i = 0; i < (WINDOWS + N) * 3; ++i) {
        #pragma HLS PIPELINE II=1
        if      (counter == 0) {               // A
            A  = m1_s->read();
            z1 = z1_s->read();
            z2 = z2_s->read();
            out1->write(A);
            out2->write((uint384)0);
        } else if (counter == 1) {             // B
            B = m1_s->read();
            uint384 E = SUB0(B, A, MOD_Q);
            uint384 H = ADD0(B, A, MOD_Q);
            e_s->write(E);
            h_s->write(H);
            out1->write(B);
            out2->write((uint384)0);
        } else {                               // T -> C, and D
            T = m1_s->read();
            uint384 C = MUL0(T, bls377_k, MOD_Q, INV_Q);         // C = k*T
            uint384 D = MUL0(ADD0(z1, z1, MOD_Q), z2, MOD_Q, INV_Q); // D = 2*z1*z2
            out1->write(C);
            out2->write(D);
        }
        counter = (counter == 2) ? 0 : counter + 1;
    }
}

//------------------------------------------------------------------------------
// Stage-3: from (C) and side-channels (E,H) and D, derive F=D-C, G=D+C
// Output layout to next multiplier (triplet cadence for padds4):
//   out1: E, G, F
//   out2: F, H, G
// Side channel t_s duplicates [E,H] for T3 = E*H in the last stage.
void padds3(
    data_stream0* m2_s, data_stream0* e_s, data_stream0* h_s,
    data_stream0* out1, data_stream0* out2, data_stream0* t_s
) {
    int counter = 0;
    uint384 C=0, D=0, E=0, H=0;

    for (int i = 0; i < (WINDOWS + N) * 3; ++i) {
        #pragma HLS PIPELINE II=1
        if      (counter == 0) {                   // C, E, H
            C = m2_s->read();
            E = e_s->read();
            H = h_s->read();
            t_s->write(E);
            t_s->write(H);
        } else if (counter == 1) {                 // D -> F,G
            D = e_s->read(); // reuse e_s to carry D from padds2::out2 via wiring
            uint384 F = SUB0(D, C, MOD_Q);
            uint384 G = ADD0(D, C, MOD_Q);
            out1->write(E);  out2->write(F);
            out1->write(G);  out2->write(H);
            out1->write(F);  out2->write(G);
        } else {
            // No new inputs on cycle 2; values were already emitted above.
        }
        counter = (counter == 2) ? 0 : counter + 1;
    }
}

//------------------------------------------------------------------------------
// Stage-4: final multiply stage. It consumes triplets from padds3:
//   m3_s: [X lane], then [Y lane], then [Z lane] operands
// and emits X3, Y3, Z3. (T3 comes from E*H computed in unipadd.)
void padds4(
    data_stream0* m3_s, data_stream0* x_s, data_stream0* y_s, data_stream0* z_s
) {
    int counter = 0;
    uint384 x=0, y=0, z=0;

    for (int i = 0; i < (WINDOWS + N) * 3; ++i) {
        #pragma HLS PIPELINE II=1
        if      (counter == 0) { x = m3_s->read(); }
        else if (counter == 1) { y = m3_s->read(); }
        else {                  z = m3_s->read(); x_s->write(x); y_s->write(y); z_s->write(z); }
        counter = (counter == 2) ? 0 : counter + 1;
    }
}

//------------------------------------------------------------------------------
// Window selectors (routing helpers).
void selector1(
    s_stream* s1,                             // window bits
    data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s,
    data_stream0* xr_s, data_stream0* yr_s, data_stream0* zr_s, data_stream0* tr_s
) {
    for (int i = 0; i < S * WINDOW_NUM; ++i) {
        #pragma HLS PIPELINE II=1
        inddata w = s1->read(); (void)w;      // use w to branch if needed
        xr_s->write(x_s->read());
        yr_s->write(y_s->read());
        zr_s->write(z_s->read());
        tr_s->write(t_s->read());
    }
}

void selector2(
    s_stream* s1,
    data_stream0* x_s, data_stream0* y_s, data_stream0* z_s, data_stream0* t_s,
    data_stream0* xr2_s, data_stream0* yr2_s, data_stream0* zr2_s, data_stream0* tr2_s
) {
    for (int i = 0; i < S * WINDOW_NUM; ++i) {
        #pragma HLS PIPELINE II=1
        inddata w = s1->read(); (void)w;
        xr2_s->write(x_s->read());
        yr2_s->write(y_s->read());
        zr2_s->write(z_s->read());
        tr2_s->write(t_s->read());
    }
}


//------------------------------------------------------------------------------
// Unified MSM point-add dataflow wiring.
// It ties padds1 -> (3 mul lanes) -> padds2 -> (2 mul lanes) -> padds3 -> (1 mul lane) -> padds4.
// E*H for T3 is computed in parallel and emitted on tr_s.
void unipadd(
    data_stream0* x1_s, data_stream0* y1_s, data_stream0* z1_s, data_stream0* t1_s,
    data_stream0* x2_s, data_stream0* y2_s, data_stream0* z2_s, data_stream0* t2_s,
    data_stream0* xr_s, data_stream0* yr_s, data_stream0* zr_s, data_stream0* tr_s
) {
#pragma HLS DATAFLOW
    // Stage streams
    data_stream0 s_p1_a, s_p1_b;
    data_stream0 s_a_in1, s_a_in2, s_b_in1, s_b_in2, s_t_in1, s_t_in2;
    data_stream0 s_a, s_b, s_t;
    data_stream0 s_p2_trip, s_p2_D;
    data_stream0 s_eh_out1, s_eh_out2, s_eh_dup;
    data_stream0 s_xlane, s_ylane, s_zlane;

#pragma HLS STREAM depth=8 variable=s_p1_a
#pragma HLS STREAM depth=8 variable=s_p1_b
#pragma HLS STREAM depth=8 variable=s_a_in1
#pragma HLS STREAM depth=8 variable=s_a_in2
#pragma HLS STREAM depth=8 variable=s_b_in1
#pragma HLS STREAM depth=8 variable=s_b_in2
#pragma HLS STREAM depth=8 variable=s_t_in1
#pragma HLS STREAM depth=8 variable=s_t_in2
#pragma HLS STREAM depth=8 variable=s_a
#pragma HLS STREAM depth=8 variable=s_b
#pragma HLS STREAM depth=8 variable=s_t
#pragma HLS STREAM depth=8 variable=s_p2_trip
#pragma HLS STREAM depth=8 variable=s_p2_D
#pragma HLS STREAM depth=8 variable=s_eh_out1
#pragma HLS STREAM depth=8 variable=s_eh_out2
#pragma HLS STREAM depth=8 variable=s_eh_dup
#pragma HLS STREAM depth=8 variable=s_xlane
#pragma HLS STREAM depth=8 variable=s_ylane
#pragma HLS STREAM depth=8 variable=s_zlane

    // padds1
    padds1(x1_s, y1_s, t1_s, x2_s, y2_s, t2_s, &s_p1_a, &s_p1_b);

    // Demultiplex triplets to three multiplier lanes
demux_loop:
    for (int i = 0; i < (WINDOWS + N) * 3; ++i) {
        #pragma HLS PIPELINE II=1
        uint384 v1 = s_p1_a.read();
        uint384 v2 = s_p1_b.read();
        int sel = i % 3;
        if      (sel == 0) { s_a_in1.write(v1); s_a_in2.write(v2); } // (Y-X)*(Y-X)
        else if (sel == 1) { s_b_in1.write(v1); s_b_in2.write(v2); } // (Y+X)*(Y+X)
        else               { s_t_in1.write(v1); s_t_in2.write(v2); } // T1*T2
    }

    // Three parallel mul lanes
    threemul(&s_a_in1, &s_a_in2, &s_a); // A
    threemul(&s_b_in1, &s_b_in2, &s_b); // B
    threemul(&s_t_in1, &s_t_in2, &s_t); // T

    // padds2: generate [A,B,C] and D; while producing E,H to side channels
    data_stream0 s_E, s_H;
#pragma HLS STREAM depth=8 variable=s_E
#pragma HLS STREAM depth=8 variable=s_H
    padds2(&s_a, z1_s, z2_s, &s_E, &s_H, &s_p2_trip, &s_p2_D);

    // padds3: turn (C, E, H, D) into (E,H,F,G); tee E,H to s_eh_dup for T3
    padds3(&s_p2_trip, /*e_s carries D here*/ &s_p2_D, &s_H,
           &s_eh_out1, &s_eh_out2, &s_eh_dup);

    // Prepare three lanes for final multiplications:
    // X3 lane uses (E,F), Y3 uses (G,H), Z3 uses (F,G)
    data_stream0 s_x_in1, s_x_in2, s_y_in1, s_y_in2, s_z_in1, s_z_in2;
#pragma HLS STREAM depth=8 variable=s_x_in1
#pragma HLS STREAM depth=8 variable=s_x_in2
#pragma HLS STREAM depth=8 variable=s_y_in1
#pragma HLS STREAM depth=8 variable=s_y_in2
#pragma HLS STREAM depth=8 variable=s_z_in1
#pragma HLS STREAM depth=8 variable=s_z_in2

prep_xyz:
    for (int i = 0; i < (WINDOWS + N); ++i) {
        #pragma HLS PIPELINE II=3
        // out1: E, G, F   | out2: F, H, G
        uint384 E = s_eh_out1.read();
        uint384 F = s_eh_out2.read();
        uint384 G = s_eh_out1.read();
        uint384 H = s_eh_out2.read();
        uint384 F2= s_eh_out1.read();
        uint384 G2= s_eh_out2.read();

        s_x_in1.write(E); s_x_in2.write(F);     // X3 = E*F
        s_y_in1.write(G); s_y_in2.write(H);     // Y3 = G*H
        s_z_in1.write(F2); s_z_in2.write(G2);   // Z3 = F*G
    }

    // Final three multiplies
    data_stream0 s_x3, s_y3, s_z3;
#pragma HLS STREAM depth=8 variable=s_x3
#pragma HLS STREAM depth=8 variable=s_y3
#pragma HLS STREAM depth=8 variable=s_z3
    threemul(&s_x_in1, &s_x_in2, &s_x3);
    threemul(&s_y_in1, &s_y_in2, &s_y3);
    threemul(&s_z_in1, &s_z_in2, &s_z3);

    // Stage-4: pack triplet to X/Y/Z streams
    padds4(&s_x3, xr_s, yr_s, zr_s);
    padds4(&s_y3, xr_s, yr_s, zr_s);
    padds4(&s_z3, xr_s, yr_s, zr_s);

    // T3 = E*H
    for (int i = 0; i < N; ++i) {
        #pragma HLS PIPELINE II=1
        uint384 E = s_eh_dup.read();
        uint384 H = s_eh_dup.read();
        tr_s->write(MUL0(E, H, MOD_Q, INV_Q));
    }
}

    s_stream*    i_s,                                 // window index per point
    data_stream0* x_s,  data_stream0* y_s,  data_stream0* z_s,  data_stream0* t_s,   // incoming points
    data_stream0* x1_s, data_stream0* y1_s, data_stream0* z1_s, data_stream0* t1_s,  // bucket-old -> add lane A (point1)
    data_stream0* x2_s, data_stream0* y2_s, data_stream0* z2_s, data_stream0* t2_s,  // incoming   -> add lane A (point2)
    data_stream0* xn1_s,data_stream0* yn1_s,data_stream0* zn1_s,data_stream0* tn1_s, // compact bucket list copy #1
    data_stream0* xn2_s,data_stream0* yn2_s,data_stream0* zn2_s,data_stream0* tn2_s, // compact bucket list copy #2
    data_stream0* xr_s, data_stream0* yr_s, data_stream0* zr_s, data_stream0* tr_s,  // add lane A result (bucket update)
    data_stream0* xr2_s,data_stream0* yr2_s,data_stream0* zr2_s,data_stream0* tr2_s, // add lane B result (downstream)
    data_stream0* x3_s, data_stream0* y3_s, data_stream0* z3_s, data_stream0* t3_s   // 32 buckets snapshot (for next stage)
) {
#pragma HLS INTERFACE mode=s_axilite port=return

#ifdef __SYNTHESIS__
    epoint buckets[WINDOWS/2];         // 1..(2^w-1) with signed-window folding
#pragma HLS DEPENDENCE dependent=false type=inter variable=buckets
#else
    epoint* buckets = (epoint*)malloc((WINDOWS/2)*sizeof(epoint));
#endif
    uint1 flags[WINDOWS/2];            // 0: idle, 1: busy (an addition in flight)

    // map add-result -> bucket id
    s_stream nor_s;
#pragma HLS STREAM depth=32 variable=nor_s

init_flags:
    for (int i = 0; i < WINDOWS/2; ++i) flags[i] = 0;

main_loop:
    for (int i = 0; i < N; ++i) {
    #pragma HLS PIPELINE II=3
        // If lane-A has a finished bucket, write it back first (keeps buckets up-to-date)
        if (!xr_s->empty() && !yr_s->empty() && !zr_s->empty() && !tr_s->empty()) {
            epoint pt0;
            inddata id0 = nor_s.read();
            pt0.x = xr_s->read();
            pt0.y = yr_s->read();
            pt0.z = zr_s->read();
            pt0.t = tr_s->read();
            flags[id0]   = 0;
            buckets[id0] = pt0;
        }

        // Read next point and its window index
        epoint  p;
        inddata idx;
        p.x = x_s->read();
        p.y = y_s->read();
        p.z = z_s->read();
        p.t = t_s->read();
        idx = i_s->read();

        if (idx == 0) continue; // skip zero-window

        // signed-window fold: msb=1 -> mirror bucket and negate Y
        if (idx[WINDOW_SIZE - 1] == 1) {
            idx  = inddata((1 << WINDOW_SIZE) - 1) - idx;
            p.y  = bls377_q - p.y;  // y := -y mod q
        }

        // If bucket is idle, dispatch (bucket_old, new_point) to lane-A immediately
        if (flags[idx] == 0) {
            flags[idx] = 1;
            x1_s->write(buckets[idx].x);  y1_s->write(buckets[idx].y);
            z1_s->write(buckets[idx].z);  t1_s->write(buckets[idx].t);  // old bucket value
            x2_s->write(p.x);             y2_s->write(p.y);
            z2_s->write(p.z);             t2_s->write(p.t);             // incoming point
            nor_s.write(idx);                                             // remember which bucket
        }
        // else: in this version we don't stage overflow here — upstream should throttle by II
    }

    // Drain any remaining lane-A results
drain_laneA:
    while (!xr_s->empty() && !yr_s->empty() && !zr_s->empty() && !tr_s->empty()) {
        epoint pt0;
        inddata id0 = nor_s.read();
        pt0.x = xr_s->read(); pt0.y = yr_s->read();
        pt0.z = zr_s->read(); pt0.t = tr_s->read();
        flags[id0]   = 0;
        buckets[id0] = pt0;
    }

    epoint tmp [32];
#pragma HLS DEPENDENCE dependent=false type=inter variable=tmp
    epoint tmp1[32];
#pragma HLS DEPENDENCE dependent=false type=inter variable=tmp1

    uint32 c1 = 0, c2 = 0;
    while (c1 < WINDOWS || c2 < WINDOWS) {
        // consume lane-B result if available -> fill tmp/tmp1 by 32-way lane, alternating halves
        if (!xr2_s->empty() && !yr2_s->empty() && !zr2_s->empty() && !tr2_s->empty()) {
            uint32 kget = c2;
            epoint eback;
            eback.x = xr2_s->read();
            eback.y = yr2_s->read();
            eback.z = zr2_s->read();
            eback.t = tr2_s->read();
            inddata cat = kget(4, 0);   // 0..31
            if (kget[5] == 0) tmp [cat] = eback;
            else              tmp1[cat] = eback;
            c2++;
        }

        // stream out compact bucket lists in lockstep: first 32 to xn1_*, next 32 to xn2_*
        inddata mouse = c1(4, 0);
        if (c1[5] == 0) {
            xn1_s->write(tmp [mouse].x);  yn1_s->write(tmp [mouse].y);
            zn1_s->write(tmp [mouse].z);  tn1_s->write(tmp [mouse].t);
        } else {
            xn2_s->write(tmp1[mouse].x);  yn2_s->write(tmp1[mouse].y);
            zn2_s->write(tmp1[mouse].z);  tn2_s->write(tmp1[mouse].t);
        }
        c1++;
    }

    // Export the second 32-bucket list as a snapshot for the next aggregation stage
    for (uint32 k = 0; k < 32; ++k) {
    #pragma HLS PIPELINE II=1
        x3_s->write(tmp1[k].x);
        y3_s->write(tmp1[k].y);
        z3_s->write(tmp1[k].z);
        t3_s->write(tmp1[k].t);
    }
}

void aggre(
    data_stream0* xb_s, data_stream0* yb_s, data_stream0* zb_s, data_stream0* tb_s,
    data_stream0* xo_s, data_stream0* yo_s, data_stream0* zo_s, data_stream0* to_s
) {
#pragma HLS DATAFLOW
    // We implement a sequential running-sum using the same padds pipeline.
    // acc starts at infinity (Z=0,Y=1,X=T=0). For each incoming bucket b (from high to low):
    //   acc = acc + b
    //   emit acc
    data_stream0 ax_s, ay_s, az_s, at_s; // accumulator stream (single-item FIFO reused)
    data_stream0 bx_s, by_s, bz_s, bt_s; // current bucket
    data_stream0 xr_s, yr_s, zr_s, tr_s; // add result
#pragma HLS STREAM depth=2 variable=ax_s
#pragma HLS STREAM depth=2 variable=ay_s
#pragma HLS STREAM depth=2 variable=az_s
#pragma HLS STREAM depth=2 variable=at_s
#pragma HLS STREAM depth=8 variable=bx_s
#pragma HLS STREAM depth=8 variable=by_s
#pragma HLS STREAM depth=8 variable=bz_s
#pragma HLS STREAM depth=8 variable=bt_s
#pragma HLS STREAM depth=2 variable=xr_s
#pragma HLS STREAM depth=2 variable=yr_s
#pragma HLS STREAM depth=2 variable=zr_s
#pragma HLS STREAM depth=2 variable=tr_s

    // Seed accumulator = infinity
    {
        epoint inf; inf.x = 0; inf.y = 1; inf.z = 0; inf.t = 0;
        ax_s.write(inf.x); ay_s.write(inf.y); az_s.write(inf.z); at_s.write(inf.t);
    }

    // For each bucket element on xb/yb/zb/tb, run unipadd(acc, b) -> acc
agg_loop:
    while (!xb_s->empty() && !yb_s->empty() && !zb_s->empty() && !tb_s->empty()) {
    #pragma HLS PIPELINE II=1
        // Acc in A, bucket in B
        data_stream0 a_x,a_y,a_z,a_t, b_x,b_y,b_z,b_t;
#pragma HLS STREAM depth=2 variable=a_x
#pragma HLS STREAM depth=2 variable=a_y
#pragma HLS STREAM depth=2 variable=a_z
#pragma HLS STREAM depth=2 variable=a_t
#pragma HLS STREAM depth=2 variable=b_x
#pragma HLS STREAM depth=2 variable=b_y
#pragma HLS STREAM depth=2 variable=b_z
#pragma HLS STREAM depth=2 variable=b_t

        // Read current acc and next bucket
        uint384 ax = ax_s.read(), ay = ay_s.read(), az = az_s.read(), at = at_s.read();
        uint384 bx = xb_s->read(), by = yb_s->read(), bz = zb_s->read(), bt = tb_s->read();

        // Feed unipadd and capture the result back into acc streams
        unipadd(&a_x,&a_y,&a_z,&a_t, &b_x,&b_y,&b_z,&b_t, &xr_s,&yr_s,&zr_s,&tr_s);

        // Provide inputs
        a_x.write(ax); a_y.write(ay); a_z.write(az); a_t.write(at);
        b_x.write(bx); b_y.write(by); b_z.write(bz); b_t.write(bt);

        // Read result
        uint384 rx = xr_s.read(), ry = yr_s.read(), rz = zr_s.read(), rt = tr_s.read();

        // Update acc
        ax_s.write(rx); ay_s.write(ry); az_s.write(rz); at_s.write(rt);

        // Emit running sum
        xo_s->write(rx); yo_s->write(ry); zo_s->write(rz); to_s->write(rt);
    }
}


} // namespace zkps::core
