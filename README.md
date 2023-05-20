# ConvertYUVtoRGB
AVS plugin for single pass convert

AVX2 (and AVX512 in future) example of single-pass YUV to RGB32 decoder (currently only point-resize for UV planes so quality is medium to draft preview).

Run several times faster at single thread in compare with AVS+ 3.7.3-times ConvertToRGB32() (looks like UV upscale + dematrix in 2 separate passes).

Currently only 4:2:0 planar input. YV12 (YUV420P8), YUV420P10, YUV420P12, YUV420P14, YUV420P16 formats (all truncated to 8bit YV12 before processing). Processing intermediates only 16bit so output precision is less in compare with ConvertToRGB32 (using 32bit intermediates). 

Usage:

DecodeYUVtoRGB(matrix=0, threads=1, gain=64, offset=0)

Default conversion do not change levels mapping (narrow YUV to narrow RGB). For 'Full' levels mapping:

DecodeYUVtoRGB(matrix=0, threads=1, gain=74, offset=-16)

Params:
1. matrix (integer), 0 - Rec.601, 1 - Rec.709, 2 - Rec.2020 (ncl ?), default 0.
2. threads (integer) - number of threads for internal OpenMP single frame multithreading, default 1.
3. offset (signed short) - value to add to 8bit RGB value after decoding, default 0 (no offset).
4. gain (short) - scaled to 64 multiplier (64 mean 1.0float - no gain), default 64 (no gain).
5. cl (bool) - cache load (true/false), default true.
6. cs (bool) - cache store (true/false), default true.
