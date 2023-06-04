# ConvertYUVtoRGB
AVS plugin for single pass convert

AVX2 (and AVX512 in future) example of single-pass YUV to RGB32 decoder (currently only point-resize for UV planes so quality is medium to draft preview).

Run several times faster at single thread in compare with AVS+ 3.7.3-times ConvertToRGB32() (looks like UV upscale + dematrix in 2 separate passes).

Currently only 4:2:0 planar input. YV12 (YUV420P8), YUV420P10, YUV420P12, YUV420P14, YUV420P16 formats (all truncated to 8bit YV12 before processing). Processing intermediates only 16bit so output precision is a bit less in compare with ConvertToRGB32 (using 32bit intermediates), also Full range may be a bit less precision in compare with narrow. 

Usage:

DecodeYUVtoRGB(matrix=0, threads=1, gain=64, offset=16) // narrow range output

Default conversion (gain 64 and offset 16) process in signed integers with zero black. After dematrix to RGB optional proc amp possible with variable gain and offset (to convert to Full/narrow/other levels mapping). 

After version 0.4.1 added 'input proc-amp' adjustments to map input range to processing domain (zero black and 219 nominal white). Defaults are for standard/narrow YUV - Ybias=-16 and UVbias=-128 and no additional scale/gain for UVs (UVgain=1.0). If input 'full-YUV' with zero black already - it is recommended to put Ybias=0 and offset=0 to keep processing and output with 'full' (if also 'full-RGB' required at output). So plugin can also do levels remapping of inputs and outputs. Inputs need to be remapped into processing internal domain of signed 16bit integer with zero black.

For 'Full' levels mapping:
DecodeYUVtoRGB(matrix=0, threads=1, gain=74, offset=0)

For non-clippig superwhites monitoring or other processing with zero black:
DecodeYUVtoRGB(matrix=0, threads=1, gain=69, offset=0)

Params:
1. matrix (integer), 0 - Rec.601, 1 - Rec.709, 2 - Rec.2020 (ncl ?), default 0.
2. threads (integer) - number of threads for internal OpenMP single frame multithreading, default 1.
3. offset (signed short) - value to add to 8bit RGB value after decoding, default 16 (offset to narrow range).
4. gain (short) - scaled to 64 multiplier (64 mean 1.0float - no gain), default 64 (no gain).
5. cl (bool) - cache load (true/false), default true.
6. cs (bool) - cache store (true/false), default false.
7. ib (integer) - 16 or 32 only allowed. 16 or 32 bit intermediate processing (also gain must be scaled to 8192 as 1.0float in 32bit intermediate), defult 16 (faster).
8. Ybias (integer) - value to add to input Y (negative mean subtraction), default -16 for narrow YUV.
9. UVbias (integer) - value to add to input U and V (negative mean subtraction), default -128 for narrow YUV (and all Digital YUV ?).
10. UVgain (float) - value to additionally scale all matrix coefficients, default 1.0f (no scale). May be used for staturation control.
