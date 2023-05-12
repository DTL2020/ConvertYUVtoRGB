# ConvertYV12toRGB
AVS plugin for single pass convert

AVX2 (and AVX512 in future) example of single-pass YV12 to RGB32 decoder (currently only point-resize for UV planes so quality is medium to draft preview).

Run several times faster at single thread in compare with AVS+ 3.7.3-times ConvertToRGB32() (looks like UV upscale + dematrix in 2 separate passes).
