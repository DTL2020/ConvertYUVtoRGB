# ConvertYV12toRGB
AVS plugin for single pass convert

AVX2 (and AVX512 in future) example of single-pass YV12 to RGB32 decoder (currently only point-resize for UV planes so quality is medium to draft preview).

Run several times faster at single thread in compare with AVS+ 3.7.3-times ConvertToRGB32() (looks like UV upscale + dematrix in 2 separate passes).

Usage:
DecodeYV12toRGB(matrix=0, threads=1)

Params:
1. matrix (integer), 0 - Rec.601, 1 - Rec.709, 2 - Rec.2020 (ncl ?), default 0.
2. threads (integer) - number of threads for internal OpenMP single frame multithreading, default 1.
