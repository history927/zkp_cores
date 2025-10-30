#include "zkps/core/ntt.hpp"
#include "zkps/config.hpp"
#include "zkps/constants.hpp"
#include "zkps/math/field.hpp"

namespace zkps::core {

using namespace zkps;
using namespace zkps::consts;
using namespace zkps::math;

// Stage-1: butterfly + twiddle multiply
void ntt1(
    data_stream* left0,  data_stream* right0,
    data_stream* left1,  data_stream* left2,
    data_stream* right1, data_stream* right2,
    data_stream* w_s,
    data_stream* outa_s, data_stream* outb_s
) {
    // First layer
    for (uint32 j = 0; j < R/2; ++j) {
        #pragma HLS PIPELINE
        uint256 ina  = left0->read();
        uint256 inb  = right0->read();
        uint256 w    = w_s->read();

        uint256 outa = ADD(ina, inb, MOD_P);             // a + b (mod p)
        uint256 outb = SUB(ina, inb, MOD_P);             // a - b (mod p)
        outb         = MUL(outb, w, MOD_P, INV_P);       // (a-b)*w (Montgomery)

        outa_s->write(outa);
        outb_s->write(outb);
    }

    // Remaining layers
    for (int i = 1; i < LOGR; ++i) {
        #pragma HLS PIPELINE off
        for (uint32 j = 0; j < R/2; ++j) {
            #pragma HLS PIPELINE
            // Select bank 1/2 based on routing bit (match your split scheme)
            bool pick2 = j & 1;
            uint256 a  = pick2 ? left2->read()  : left1->read();
            uint256 b  = pick2 ? right2->read() : right1->read();
            uint256 w  = w_s->read();

            uint256 outa = ADD(a, b, MOD_P);
            uint256 outb = SUB(a, b, MOD_P);
            outb         = MUL(outb, w, MOD_P, INV_P);

            outa_s->write(outa);
            outb_s->write(outb);
        }
    }
}

// Stage-2: route outa/outb back to left/right banks and final pass-through
void ntt2(
    data_stream* outa_s, data_stream* outb_s,
    data_stream* left1,  data_stream* left2,
    data_stream* right1, data_stream* right2,
    data_stream* c1,     data_stream* c2
) {
    // Route back for the first LOGR-1 layers
    for (int i = 0; i < LOGR - 1; ++i) {
        #pragma HLS PIPELINE off
        for (uint32 j = 0; j < R/2; ++j) {
            #pragma HLS PIPELINE
            uint256 a = outa_s->read();
            uint256 b = outb_s->read();
            if (((j >> i) & 1) == 0) {
                left1->write(a);
                left2->write(b);
            } else {
                right1->write(a);
                right2->write(b);
            }
        }
    }

    // Final pass-through to contiguous outputs
    for (uint32 k = 0; k < R/2; ++k) {
        c1->write(outa_s->read());
        c2->write(outb_s->read());
    }
}

} // namespace zkps::core
