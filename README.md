# ZKPS (HLS) — Vitis-ready

Self-contained HLS code for NTT and MSM streaming pipelines.  
The files under **`vitis/`** (`zkps.cpp`, `zkps.hpp`) are the **merged, single translation unit** intended for **`v++`**. All other sources are split for development, readability, and unit testing.

---

## Project layout

```

include/zkps/
config.hpp           # sizes: R, LOGR, S, WINDOW_*
constants.hpp        # field/curve constants & Montgomery inverses
types.hpp            # ap_uint typedefs, streams, epoint, etc.
core/
    msm.hpp            # MSM pipeline API (padds/selector/bucket)
    ntt.hpp            # NTT stage APIs (ntt1 / ntt2)
io/
    io.hpp             # read_*/writeback/lazy_* loaders & small cores
math/
    field.hpp          # mod add/sub, Montgomery mul (256/384)
    karatsuba.hpp      # 128/256/384 multiplies (divide & conquer)
    montgomery.hpp     # 256/384 Montgomery reduction


src/
    core/
        msm.cpp            # threemul, padds1..4, selector1/2, bucket_part, unipadd
        ntt.cpp            # ntt1 / ntt2
    io/
        io.cpp             # read_*, writeback, subcore/mulcore, lazy_*, data_loader_*
    math/
        math.cpp           # field.hpp implementation (ADD/SUB/MUL wrappers)
        montgomery.cpp     # montgomery_reduce_256/384

vitis/
    zkps.cpp             # merged single-TU source for v++
    zkps.hpp             # matching header for the merged build

````



## Vitis / `v++` flow (use `vitis/`)

Compile kernels from the **merged** sources:

```bash
# Compile to XO
v++ -c -t hw \
  --platform <platform> \
  -k <kernel_name> \
  -I include \
  vitis/zkps.cpp \
  -o <kernel_name>.xo
```

```bash
# Link to XCLBIN
v++ -l -t hw \
  --platform <platform> \
  <kernel_name>.xo \
  -o app.xclbin
```



