#include "include\avisynth.h"
#include <windows.h>
#include <immintrin.h>

#pragma warning(disable : 4309)

#define RGB_DIV_SHIFT 6
#define RGB_DIV_SHIFT_32 13

#define RND_16 32 // 1<<3 ?
#define RND_32 4096 // 1<<12 ?

typedef enum ColorRange_e {
	AVS_RANGE_FULL = 0,
	AVS_RANGE_LIMITED = 1
} ColorRange_e;

// https://www.itu.int/rec/T-REC-H.265-202108-I
/* ITU-T H.265 Table E.5 */
typedef enum Matrix_e {
	AVS_MATRIX_RGB = 0, /* The identity matrix. Typically used for RGB, may also be used for XYZ */
	AVS_MATRIX_BT709 = 1, /* ITU-R Rec. BT.709-5 */
	AVS_MATRIX_UNSPECIFIED = 2, /* Image characteristics are unknown or are determined by the application */
	AVS_MATRIX_BT470_M = 4, // instead of AVS_MATRIX_FCC
	// FCC Title 47 Code of Federal Regulations (2003) 73.682 (a) (20)
	// Rec. ITU-R BT.470-6 System M (historical)
	AVS_MATRIX_BT470_BG = 5, /* Equivalent to 6. */
	// ITU-R Rec. BT.470-6 System B, G (historical)
	// Rec. ITU-R BT.601-7 625
	// Rec. ITU-R BT.1358-0 625 (historical)
	// Rec. ITU-R BT.1700-0 625 PAL and 625 SECAM
	AVS_MATRIX_ST170_M = 6,  /* Equivalent to 5. */
	// Rec. ITU-R BT.601-7 525
	// Rec. ITU-R BT.1358-1 525 or 625 (historical)
	// Rec. ITU-R BT.1700-0 NTSC
	// SMPTE ST 170 (2004)
	// SMPTE 170M (2004)
	AVS_MATRIX_ST240_M = 7, // SMPTE ST 240 (1999, historical)
	AVS_MATRIX_YCGCO = 8,
	AVS_MATRIX_BT2020_NCL = 9,
	// Rec. ITU-R BT.2020 non-constant luminance system
	// Rec. ITU-R BT.2100-2 Y'CbCr
	AVS_MATRIX_BT2020_CL = 10, /* Rec. ITU-R BT.2020 constant luminance system */
	AVS_MATRIX_CHROMATICITY_DERIVED_NCL = 12, /* Chromaticity derived non-constant luminance system */
	AVS_MATRIX_CHROMATICITY_DERIVED_CL = 13, /* Chromaticity derived constant luminance system */
	AVS_MATRIX_ICTCP = 14, // REC_2100_ICTCP, Rec. ITU-R BT.2100-2 ICTCP
	AVS_MATRIX_AVERAGE = 9999, // Avisynth compatibility
} Matrix_e;


struct ConversionMatrix {
	int y_r, y_g, y_b;
	// for grayscale conversion these may not needed
	int u_r, u_g, u_b;
	int v_r, v_g, v_b;

	// used in YUY2 RGB->YUY2 asm
	int ku, ku_luma;
	int kv, kv_luma;

	float y_r_f, y_g_f, y_b_f;
	float u_r_f, u_g_f, u_b_f;
	float v_r_f, v_g_f, v_b_f;

	int offset_y;
	float offset_y_f;
};

// 8 bit fullscale to float
static AVS_FORCEINLINE float c8tof(int color) {
	return color / 255.0f;
}

static AVS_FORCEINLINE float uv8tof(int color) {
	return (color - 128) / 255.f; // consistent with convert_uintN_to_float_c
}

static void BuildMatrix_Yuv2Rgb_core(double Kr, double Kb, int int_arith_shift, bool full_scale, int bits_per_pixel, ConversionMatrix& matrix)
{
	int Sy, Suv, Oy;
	float Sy_f, Suv_f, Oy_f;

	if (bits_per_pixel <= 16) {
		Oy = full_scale ? 0 : (16 << (bits_per_pixel - 8));
		Oy_f = (float)Oy; // for 16 bits

		int ymin = (full_scale ? 0 : 16) << (bits_per_pixel - 8);
		int max_pixel_value = (1 << bits_per_pixel) - 1;
		int ymax = full_scale ? max_pixel_value : (235 << (bits_per_pixel - 8));
		Sy = ymax - ymin;
		Sy_f = (float)Sy;

		int cmin = full_scale ? 0 : (16 << (bits_per_pixel - 8));
		int cmax = full_scale ? max_pixel_value : (240 << (bits_per_pixel - 8));
		Suv = (cmax - cmin) / 2;
		Suv_f = (cmax - cmin) / 2.0f;
	}
	else {
		Oy_f = full_scale ? 0.0f : (16.0f / 255.0f);
		Oy = full_scale ? 0 : 16; // n/a

		Sy_f = full_scale ? c8tof(255) : (c8tof(235) - c8tof(16));
		Suv_f = full_scale ? (0.5f - -0.5f) / 2 : (uv8tof(240) - uv8tof(16)) / 2;
	}


	/*
	  Kr   = {0.299, 0.2126}
	  Kb   = {0.114, 0.0722}
	  Kg   = 1 - Kr - Kb // {0.587, 0.7152}
	  Srgb = 255
	  Sy   = {219, 255}   // { 235-16, 255-0 }
	  Suv  = {112, 127}   // { (240-16)/2, (255-0)/2 }
	  Oy   = {16, 0}
	  Ouv  = 128

	  Y =(y-Oy)  / Sy                         // 0..1
	  U =(u-Ouv) / Suv                        //-1..1
	  V =(v-Ouv) / Suv

	  R = Y                  + V*(1-Kr)       // 0..1
	  G = Y - U*(1-Kb)*Kb/Kg - V*(1-Kr)*Kr/Kg
	  B = Y + U*(1-Kb)

	  r = R*Srgb                              // 0..255   0..65535
	  g = G*Srgb
	  b = B*Srgb
	*/

	const double mulfac = double(1 << int_arith_shift); // integer aritmetic precision scale

	const double Kg = 1. - Kr - Kb;

	if (bits_per_pixel <= 16) {
		const int Srgb = (1 << bits_per_pixel) - 1;  // 255;
		matrix.y_b = (int)(Srgb * 1.000 * mulfac / Sy + 0.5); //Y
		matrix.u_b = (int)(Srgb * (1 - Kb) * mulfac / Suv + 0.5); //U
		matrix.v_b = (int)(Srgb * 0.000 * mulfac / Suv + 0.5); //V
		matrix.y_g = (int)(Srgb * 1.000 * mulfac / Sy + 0.5);
		matrix.u_g = (int)(Srgb * (Kb - 1) * Kb / Kg * mulfac / Suv + 0.5);
		matrix.v_g = (int)(Srgb * (Kr - 1) * Kr / Kg * mulfac / Suv + 0.5);
		matrix.y_r = (int)(Srgb * 1.000 * mulfac / Sy + 0.5);
		matrix.u_r = (int)(Srgb * 0.000 * mulfac / Suv + 0.5);
		matrix.v_r = (int)(Srgb * (1 - Kr) * mulfac / Suv + 0.5);
		matrix.offset_y = -Oy;
	}

	double Srgb_f = bits_per_pixel == 32 ? 1.0 : ((1 << bits_per_pixel) - 1);
	matrix.y_b_f = (float)(Srgb_f * 1.000 / Sy_f); //Y
	matrix.u_b_f = (float)(Srgb_f * (1 - Kb) / Suv_f); //U
	matrix.v_b_f = (float)(Srgb_f * 0.000 / Suv_f); //V
	matrix.y_g_f = (float)(Srgb_f * 1.000 / Sy_f);
	matrix.u_g_f = (float)(Srgb_f * (Kb - 1) * Kb / Kg / Suv_f);
	matrix.v_g_f = (float)(Srgb_f * (Kr - 1) * Kr / Kg / Suv_f);
	matrix.y_r_f = (float)(Srgb_f * 1.000 / Sy_f);
	matrix.u_r_f = (float)(Srgb_f * 0.000 / Suv_f);
	matrix.v_r_f = (float)(Srgb_f * (1 - Kr) / Suv_f);
	matrix.offset_y_f = -Oy_f;
}

bool do_BuildMatrix_Yuv2Rgb(int _Matrix, int _ColorRange, int int_arith_shift, int bits_per_pixel, ConversionMatrix& matrix)
{
	if (_ColorRange != ColorRange_e::AVS_RANGE_FULL && _ColorRange != ColorRange_e::AVS_RANGE_LIMITED)
		return false;
	const bool is_full = _ColorRange == ColorRange_e::AVS_RANGE_FULL;

	if (_Matrix == Matrix_e::AVS_MATRIX_BT470_BG || _Matrix == Matrix_e::AVS_MATRIX_ST170_M) { // 601
		BuildMatrix_Yuv2Rgb_core(0.299,  /* 0.587  */ 0.114, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_BT709) {
		BuildMatrix_Yuv2Rgb_core(0.2126, /* 0.7152 */ 0.0722, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_AVERAGE) { // non-standard!
		BuildMatrix_Yuv2Rgb_core(1.0 / 3, /* 1.0/3 */ 1.0 / 3, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_BT2020_CL || _Matrix == Matrix_e::AVS_MATRIX_BT2020_NCL) {
		BuildMatrix_Yuv2Rgb_core(0.2627, /* 0.6780 */ 0.0593, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_BT470_M) {
		BuildMatrix_Yuv2Rgb_core(0.3, /* 0.59 */ 0.11, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_ST240_M) {
		BuildMatrix_Yuv2Rgb_core(0.212, /* 0.701 */ 0.087, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_RGB) {
		BuildMatrix_Yuv2Rgb_core(0.0, /*  */ 0.0, int_arith_shift, is_full, bits_per_pixel, matrix);
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_ICTCP) {
		// not supported REC_2100_LMS
		return false;
	}
	else if (_Matrix == Matrix_e::AVS_MATRIX_YCGCO) {
		// not supported
		return false;
	}
	else {
		return false;
	}
	return true;
}

class DecodeYUVtoRGB : public GenericVideoFilter
{
	int threads;
	int _cpuFlags;

	short Kr; // div64 for 16bit imm, div8192 for 32bit imm
	short Kb; // div64 for 16bit imm, div8192 for 32bit imm
	short Kgu; // div64 for 16bit imm, div8192 for 32bit imm
	short Kgv; // div64 for 16bit imm, div8192 for 32bit imm

	int Matrix;

	short RGBgain;
	short RGBoffset;

	bool bCacheLoad;
	bool bCacheStore;

	int iImmBits;

public:
	DecodeYUVtoRGB(PClip _child, int threads_, int _matrix, int _gain, int _offset, bool _cl, bool _cs, int _ImmBits, IScriptEnvironment* env) : GenericVideoFilter(_child),
		threads(threads_),
		Matrix(_matrix),
		RGBgain((short)_gain),
		RGBoffset((short)_offset),
		bCacheLoad(_cl),
		bCacheStore(_cs),
		iImmBits(_ImmBits)
	{
		_cpuFlags = env->GetCPUFlags();

//		vi.pixel_type = VideoInfo::CS_RGBP8; // temp perf test
		vi.pixel_type = VideoInfo::CS_BGR32;

		if (!(_cpuFlags & CPUF_AVX2))
		{
			env->ThrowError("DecodeYUVtoRGB: Only AVX2 and later SIMD co-processor supported.");
		}

		if ((iImmBits != 16) && (iImmBits != 32))
		{
			env->ThrowError("DecodeYUVtoRGB: Intermediate bits precision may be only 16 or 32.");
		}

		if (iImmBits == 16)
		{

			if (Matrix == 0) // 601
			{
				Kr = 90; // Kr of 601 ? , div64  1.402
				Kb = 113; // Kb of 601 ? , div64  1.772
				Kgu = 22; // Kgu of 601 ? , div64 0.344
				Kgv = 46; // Kgv of 601 ? , div64 0.714 */
			}
			else if (Matrix == 1) // 709
			{
				Kr = 100; // Kr of 709 ? , div64  1.575
				Kb = 119; // Kb of 709 ? , div64  1.856
				Kgu = 12; // Kgu of 709 ? , div64 0.187
				Kgv = 30; // Kgv of 709 ? , div64 0.468
			}
			else if (Matrix == 2) // 2020 (ncl ?)
			{
				Kr = 94; // Kr of 2020 ? , div64  1.475
				Kb = 120; // Kb of 2020 ? , div64  1.88
				Kgu = 10; // Kgu of 2020 ? , div64 0.165
				Kgv = 37; // Kgv of 2020 ? , div64 0.571
			}
			else
				env->ThrowError("DecodeYUVtoRGB: matrix %d not supported.", Matrix);
		}

		if (iImmBits == 32)
		{

			if (Matrix == 0) // 601
			{
				Kr = 11485; // Kr of 601 ? , div8192  1.402
				Kb = 14516; // Kb of 601 ? , div8192  1.772
				Kgu = 2819; // Kgu of 601 ? , div8192 0.344136
				Kgv = 5850; // Kgv of 601 ? , div8192 0.714136 */
			}
			else if (Matrix == 1) // 709
			{
				Kr = 12901; // Kr of 709 ? , div8192  1.5748
				Kb = 15201; // Kb of 709 ? , div8192  1.8556
				Kgu = 1535; // Kgu of 709 ? , div8192 0.187324
				Kgv = 3835; // Kgv of 709 ? , div8192 0.468124
			}
			else if (Matrix == 2) // 2020 (ncl ?)
			{
				Kr = 12081; // Kr of 2020 ? , div8192  1.4747
				Kb = 15412; // Kb of 2020 ? , div8192  1.8814
				Kgu = 1348; // Kgu of 2020 ? , div8192 0.164553
				Kgv = 4681; // Kgv of 2020 ? , div8192 0.571392
			}
			else
				env->ThrowError("DecodeYUVtoRGB: matrix %d not supported.", Matrix);
		}

	}

	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env)
	{
		PVideoFrame dst = env->NewVideoFrame(vi);
		VideoInfo vi_src = child->GetVideoInfo();
		PVideoFrame src = child->GetFrame(n, env);

		DWORD dwOldProtY;
		DWORD dwOldProtU;
		DWORD dwOldProtV;
		bool bRes;

		if (iImmBits == 16)
		{

			if (vi_src.ComponentSize() == 1)
			{

				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm16<true, true, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm16<false, true, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<true, false, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<false, false, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 10))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm16<true, true, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm16<false, true, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<true, false, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<false, false, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 12))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm16<true, true, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm16<false, true, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<true, false, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<false, false, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 14))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm16<true, true, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm16<false, true, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<true, false, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<false, false, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 16))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm16<true, true, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm16<false, true, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<true, false, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm16<false, false, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else
				env->ThrowError("DecodeYV12toRGB: Only 8bit to 16bit input supported.");
		}

		if (iImmBits == 32)
		{
			if (vi_src.ComponentSize() == 1)
			{

				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm32<true, true, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm32<false, true, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<true, false, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<false, false, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			}
/*			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 10))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm32<true, true, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm32<false, true, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<true, false, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<false, false, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 12))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm32<true, true, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm32<false, true, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<true, false, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<false, false, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 14))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm32<true, true, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm32<false, true, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<true, false, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<false, false, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}
			else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 16))
			{
				if (bCacheLoad && bCacheStore)
					DecodeYUV420imm32<true, true, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && bCacheStore)
					DecodeYUV420imm32<false, true, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<true, false, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

				if (!bCacheLoad && !bCacheStore)
					DecodeYUV420imm32<false, false, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
			}*/
			else
				env->ThrowError("DecodeYV12toRGB: Only 8bit to 16bit input supported.");
		}
		return dst;
	}

	template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
	void DecodeYUV420imm16(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags);

	template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
	void DecodeYUV420imm32(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags);

	template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
	void DecodeYUV420imm32_mat(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags);

};


template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
void DecodeYUVtoRGB::DecodeYUV420imm16(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags)
{
	auto srcp_Y = src->GetReadPtr(PLANAR_Y);
	auto srcp_U = src->GetReadPtr(PLANAR_U);
	auto srcp_V = src->GetReadPtr(PLANAR_V);


	auto height = src->GetHeight(PLANAR_Y);
	auto row_size_Y = src->GetRowSize(PLANAR_Y);
	auto row_size_U = src->GetRowSize(PLANAR_U);
	auto row_size_V = src->GetRowSize(PLANAR_V);
	auto src_pitch_Y = src->GetPitch(PLANAR_Y);
	auto src_pitch_U = src->GetPitch(PLANAR_U);
	auto src_pitch_V = src->GetPitch(PLANAR_V);

	auto dstp_R = dst->GetWritePtr(4);
	auto dstp_G = dst->GetWritePtr(6);
	auto dstp_B = dst->GetWritePtr(2);

	auto dstp_BGRA = dst->GetWritePtr(2);
	auto dst_pitch_BGRA = dst->GetPitch();

	auto dst_pitch_R = dst->GetPitch(PLANAR_R);
	auto dst_pitch_G = dst->GetPitch(PLANAR_G);
	auto dst_pitch_B = dst->GetPitch(PLANAR_B);


#pragma omp parallel for num_threads(threads)
	for (uint64_t y = 0; y < height; y++)
	{
		unsigned char* l_dstp_R = dstp_R + y * dst_pitch_R;
		unsigned char* l_dstp_G = dstp_G + y * dst_pitch_G;
		unsigned char* l_dstp_B = dstp_B + y * dst_pitch_B;

		unsigned char* l_dstp_BGRA = dstp_BGRA + (height - y - 1) * dst_pitch_BGRA; // reverse scan for RGB interleaved

		unsigned char* l_srcp_Y = (unsigned char*)srcp_Y + y * src_pitch_Y;
		unsigned char* l_srcp_U = (unsigned char*)srcp_U + (y >> 1) * src_pitch_U;
		unsigned char* l_srcp_V = (unsigned char*)srcp_V + (y >> 1) * src_pitch_V;


			if (cpuFlags & CPUF_AVX2) // use AVX2
			{
				int col64;
				int row_size_Ysamples;
				if (cs == 1) // 8bit
				{
					col64 = row_size_Y - (row_size_Y % 64);
					row_size_Ysamples = row_size_Y;
				}
				else // 10 to 16 bit
				{
					col64 = (row_size_Y / 2) - ((row_size_Y / 2) % 64);
					row_size_Ysamples = row_size_Y / 2;
				}


/*				int row_proc_size;
				if (cs == 1) // 8bit
					row_proc_size = row_size_Y;
				else // 10 to 16 bit
					row_proc_size = (row_size_Y / 2); 
					*/
				__m256i ymm_w_cbias = _mm256_set1_epi16(128);
				__m256i ymm_w_rounder = _mm256_set1_epi16(RND_16);

				__m256i ymm_wKr = _mm256_set1_epi16(Kr); 
				__m256i ymm_wKb = _mm256_set1_epi16(Kb); 

				__m256i ymm_wKgu = _mm256_set1_epi16(Kgu); 
				__m256i ymm_wKgv = _mm256_set1_epi16(Kgv);

				int UVpref = 0;

				
				// fill Y with 0 to 63 - debug
				/*
				for (int idx = 0; idx < 64; idx++)
				{
					l_srcp_Y[idx] = (unsigned int)idx;
				}
				*/
/*				unsigned short* lus_srcp_Y = (unsigned short*)l_srcp_Y;
				for (int idx = 0; idx < 64; idx++)
				{
					lus_srcp_Y[idx] = (unsigned short)(idx * 4);
				}
	*/			
				/*
				// fill U with 0 to 32 - debug
				for (int idx = 0; idx < 32; idx++)
				{
					l_srcp_U[idx] = (unsigned char)idx;
				}

				// fill V with //0 to 32 - debug
				for (int idx = 0; idx < 32; idx++)
				{
					l_srcp_V[idx] = 128;
				}
				*/
/*				unsigned short* lus_srcp_U = (unsigned short*)l_srcp_U;
				for (int idx = 0; idx < 32; idx++)
				{
					lus_srcp_U[idx] = (unsigned short)(idx * 4);
				}
*/				
							   
				for (int col = 0; col < col64; col += 64)
				{
					__m256i ymm_Y0_16l;
					__m256i ymm_Y1_16l;

					__m256i ymm_Y0_16h;
					__m256i ymm_Y1_16h;

					__m256i ymm_U0_16l;
					__m256i ymm_V0_16l;

					__m256i ymm_U0_16h;
					__m256i ymm_V0_16h;

					__m256i ymm_U1_16l;
					__m256i ymm_V1_16l;

					__m256i ymm_U1_16h;
					__m256i ymm_V1_16h;


					if (cs == 1)
					{

						if (!bCacheLoad)
						{
							_mm_prefetch((const CHAR*)(l_srcp_Y + 64), _MM_HINT_NTA);
							if (UVpref % 2 == 0)
							{
								_mm_prefetch((const CHAR*)(l_srcp_U + 64), _MM_HINT_NTA);
								_mm_prefetch((const CHAR*)(l_srcp_V + 64), _MM_HINT_NTA);
							}

							UVpref++;
						}

						__m256i ymm0_Y0 = _mm256_load_si256((const __m256i*)l_srcp_Y); // should always load from 64-bit aligned start of row in AVS+ 3.7.3 (and later ?)
						__m256i ymm1_Y1 = _mm256_load_si256((const __m256i*)(l_srcp_Y + 32));
						__m256i ymm2_U = _mm256_load_si256((const __m256i*)(l_srcp_U));
						__m256i ymm3_V = _mm256_load_si256((const __m256i*)(l_srcp_V));


						ymm_Y0_16l = _mm256_unpacklo_epi8(ymm0_Y0, _mm256_setzero_si256());
						ymm_Y1_16l = _mm256_unpacklo_epi8(ymm1_Y1, _mm256_setzero_si256());

						ymm_Y0_16h = _mm256_unpackhi_epi8(ymm0_Y0, _mm256_setzero_si256());
						ymm_Y1_16h = _mm256_unpackhi_epi8(ymm1_Y1, _mm256_setzero_si256());

						__m256i ymm_U_dl = _mm256_unpacklo_epi8(ymm2_U, ymm2_U);
						__m256i ymm_V_dl = _mm256_unpacklo_epi8(ymm3_V, ymm3_V);

						__m256i ymm_U_dh = _mm256_unpackhi_epi8(ymm2_U, ymm2_U);
						__m256i ymm_V_dh = _mm256_unpackhi_epi8(ymm3_V, ymm3_V);

						__m256i ymm_U0 = _mm256_permute2x128_si256(ymm_U_dl, ymm_U_dh, 0x20);
						__m256i ymm_V0 = _mm256_permute2x128_si256(ymm_V_dl, ymm_V_dh, 0x20);

						__m256i ymm_U1 = _mm256_permute2x128_si256(ymm_U_dl, ymm_U_dh, 0x31);
						__m256i ymm_V1 = _mm256_permute2x128_si256(ymm_V_dl, ymm_V_dh, 0x31);

						ymm_U0_16l = _mm256_unpacklo_epi8(ymm_U0, _mm256_setzero_si256());
						ymm_V0_16l = _mm256_unpacklo_epi8(ymm_V0, _mm256_setzero_si256());
						
						ymm_U0_16h = _mm256_unpackhi_epi8(ymm_U0, _mm256_setzero_si256());
						ymm_V0_16h = _mm256_unpackhi_epi8(ymm_V0, _mm256_setzero_si256());

						ymm_U1_16l = _mm256_unpacklo_epi8(ymm_U1, _mm256_setzero_si256());
						ymm_V1_16l = _mm256_unpacklo_epi8(ymm_V1, _mm256_setzero_si256());

						ymm_U1_16h = _mm256_unpackhi_epi8(ymm_U1, _mm256_setzero_si256());
						ymm_V1_16h = _mm256_unpackhi_epi8(ymm_V1, _mm256_setzero_si256());

					}
					if (cs == 2)
					{
						if (!bCacheLoad)
						{
							_mm_prefetch((const CHAR*)(l_srcp_Y + 128), _MM_HINT_NTA);
 							_mm_prefetch((const CHAR*)(l_srcp_Y + 128 + 64), _MM_HINT_NTA); 

							_mm_prefetch((const CHAR*)(l_srcp_U + 64), _MM_HINT_NTA);
							_mm_prefetch((const CHAR*)(l_srcp_V + 64), _MM_HINT_NTA);
						}

						__m256i ymm_Y0_l = _mm256_load_si256((const __m256i*)l_srcp_Y); // should always load from 64-bit aligned start of row in AVS+ 3.7.3 (and later ?)
						__m256i ymm_Y0_h = _mm256_load_si256((const __m256i*)(l_srcp_Y + 32));
						__m256i ymm_Y1_l = _mm256_load_si256((const __m256i*)(l_srcp_Y + 64));
						__m256i ymm_Y1_h = _mm256_load_si256((const __m256i*)(l_srcp_Y + 96));

						__m256i ymm_U_l = _mm256_load_si256((const __m256i*)(l_srcp_U));
						__m256i ymm_U_h = _mm256_load_si256((const __m256i*)(l_srcp_U + 32));
						
						__m256i ymm_V_l = _mm256_load_si256((const __m256i*)(l_srcp_V));
						__m256i ymm_V_h = _mm256_load_si256((const __m256i*)(l_srcp_V + 32));

						ymm_Y0_16l = _mm256_permute2x128_si256(ymm_Y0_l, ymm_Y0_h, 0x20);
						ymm_Y0_16h = _mm256_permute2x128_si256(ymm_Y0_l, ymm_Y0_h, 0x31);

						ymm_Y1_16l = _mm256_permute2x128_si256(ymm_Y1_l, ymm_Y1_h, 0x20);
						ymm_Y1_16h = _mm256_permute2x128_si256(ymm_Y1_l, ymm_Y1_h, 0x31);

						if (bps == 10)
						{
							ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 2);
							ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 2);

							ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 2);
							ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 2);

							ymm_U_l = _mm256_srli_epi16(ymm_U_l, 2);
							ymm_U_h = _mm256_srli_epi16(ymm_U_h, 2);

							ymm_V_l = _mm256_srli_epi16(ymm_V_l, 2);
							ymm_V_h = _mm256_srli_epi16(ymm_V_h, 2);
						}

						if (bps == 12)
						{
							ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 4);
							ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 4);

							ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 4);
							ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 4);

							ymm_U_l = _mm256_srli_epi16(ymm_U_l, 4);
							ymm_U_h = _mm256_srli_epi16(ymm_U_h, 4);

							ymm_V_l = _mm256_srli_epi16(ymm_V_l, 4);
							ymm_V_h = _mm256_srli_epi16(ymm_V_h, 4);
						}

						if (bps == 14)
						{
							ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 6);
							ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 6);

							ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 6);
							ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 6);

							ymm_U_l = _mm256_srli_epi16(ymm_U_l, 6);
							ymm_U_h = _mm256_srli_epi16(ymm_U_h, 6);

							ymm_V_l = _mm256_srli_epi16(ymm_V_l, 6);
							ymm_V_h = _mm256_srli_epi16(ymm_V_h, 6);
						}

						if (bps == 16)
						{
							ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 8);
							ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 8);

							ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 8);
							ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 8);

							ymm_U_l = _mm256_srli_epi16(ymm_U_l, 8);
							ymm_U_h = _mm256_srli_epi16(ymm_U_h, 8);

							ymm_V_l = _mm256_srli_epi16(ymm_V_l, 8);
							ymm_V_h = _mm256_srli_epi16(ymm_V_h, 8);
						}


						ymm_U0_16l = _mm256_unpacklo_epi16(ymm_U_l, ymm_U_l);
						ymm_V0_16l = _mm256_unpacklo_epi16(ymm_V_l, ymm_V_l);

						ymm_U1_16l = _mm256_unpacklo_epi16(ymm_U_h, ymm_U_h);
						ymm_V1_16l = _mm256_unpacklo_epi16(ymm_V_h, ymm_V_h);

						ymm_U0_16h = _mm256_unpackhi_epi16(ymm_U_l, ymm_U_l);
						ymm_V0_16h = _mm256_unpackhi_epi16(ymm_V_l, ymm_V_l);

						ymm_U1_16h = _mm256_unpackhi_epi16(ymm_U_h, ymm_U_h);
						ymm_V1_16h = _mm256_unpackhi_epi16(ymm_V_h, ymm_V_h);


					}

					ymm_U0_16l = _mm256_sub_epi16(ymm_U0_16l, ymm_w_cbias); // superscalarity of 3 at IceLake/SkyLake
					ymm_V0_16l = _mm256_sub_epi16(ymm_V0_16l, ymm_w_cbias);

					ymm_U0_16h = _mm256_sub_epi16(ymm_U0_16h, ymm_w_cbias);
					ymm_V0_16h = _mm256_sub_epi16(ymm_V0_16h, ymm_w_cbias);

					ymm_U1_16l = _mm256_sub_epi16(ymm_U1_16l, ymm_w_cbias);
					ymm_V1_16l = _mm256_sub_epi16(ymm_V1_16l, ymm_w_cbias);

					ymm_U1_16h = _mm256_sub_epi16(ymm_U1_16h, ymm_w_cbias);
					ymm_V1_16h = _mm256_sub_epi16(ymm_V1_16h, ymm_w_cbias);


					__m256i ymm_V_dl16l_addR = _mm256_mullo_epi16(ymm_V0_16l, ymm_wKr);
					__m256i ymm_V_dl16h_addR = _mm256_mullo_epi16(ymm_V0_16h, ymm_wKr);

					__m256i ymm_V_dh16l_addR = _mm256_mullo_epi16(ymm_V1_16l, ymm_wKr);
					__m256i ymm_V_dh16h_addR = _mm256_mullo_epi16(ymm_V1_16h, ymm_wKr);

					
					__m256i ymm_U_dl16l_addB = _mm256_mullo_epi16(ymm_U0_16l, ymm_wKb);
					__m256i ymm_U_dl16h_addB = _mm256_mullo_epi16(ymm_U0_16h, ymm_wKb);

					__m256i ymm_U_dh16l_addB = _mm256_mullo_epi16(ymm_U1_16l, ymm_wKb);
					__m256i ymm_U_dh16h_addB = _mm256_mullo_epi16(ymm_U1_16h, ymm_wKb);

					ymm_V_dl16l_addR = _mm256_add_epi16(ymm_V_dl16l_addR, ymm_w_rounder);
					ymm_V_dl16h_addR = _mm256_add_epi16(ymm_V_dl16h_addR, ymm_w_rounder);

					ymm_V_dh16l_addR = _mm256_add_epi16(ymm_V_dh16l_addR, ymm_w_rounder);
					ymm_V_dh16h_addR = _mm256_add_epi16(ymm_V_dh16h_addR, ymm_w_rounder);

					ymm_U_dl16l_addB = _mm256_add_epi16(ymm_U_dl16l_addB, ymm_w_rounder);
					ymm_U_dl16h_addB = _mm256_add_epi16(ymm_U_dl16h_addB, ymm_w_rounder);

					ymm_U_dh16l_addB = _mm256_add_epi16(ymm_U_dh16l_addB, ymm_w_rounder);
					ymm_U_dh16h_addB = _mm256_add_epi16(ymm_U_dh16h_addB, ymm_w_rounder);


					ymm_V_dl16l_addR = _mm256_srai_epi16(ymm_V_dl16l_addR, 6);
					ymm_V_dl16h_addR = _mm256_srai_epi16(ymm_V_dl16h_addR, 6);

					ymm_V_dh16l_addR = _mm256_srai_epi16(ymm_V_dh16l_addR, 6);
					ymm_V_dh16h_addR = _mm256_srai_epi16(ymm_V_dh16h_addR, 6);


					ymm_U_dl16l_addB = _mm256_srai_epi16(ymm_U_dl16l_addB, 6);
					ymm_U_dl16h_addB = _mm256_srai_epi16(ymm_U_dl16h_addB, 6);

					ymm_U_dh16l_addB = _mm256_srai_epi16(ymm_U_dh16l_addB, 6);
					ymm_U_dh16h_addB = _mm256_srai_epi16(ymm_U_dh16h_addB, 6);


					// R
					__m256i ymm_R_0_16l = _mm256_add_epi16(ymm_Y0_16l, ymm_V_dl16l_addR);
					__m256i ymm_R_0_16h = _mm256_add_epi16(ymm_Y0_16h, ymm_V_dl16h_addR);

					__m256i ymm_R_1_16l = _mm256_add_epi16(ymm_Y1_16l, ymm_V_dh16l_addR);
					__m256i ymm_R_1_16h = _mm256_add_epi16(ymm_Y1_16h, ymm_V_dh16h_addR);

					// B
					__m256i ymm_B_0_16l = _mm256_add_epi16(ymm_Y0_16l, ymm_U_dl16l_addB);
					__m256i ymm_B_0_16h = _mm256_add_epi16(ymm_Y0_16h, ymm_U_dl16h_addB);

					__m256i ymm_B_1_16l = _mm256_add_epi16(ymm_Y1_16l, ymm_U_dh16l_addB);
					__m256i ymm_B_1_16h = _mm256_add_epi16(ymm_Y1_16h, ymm_U_dh16h_addB);

					//G
					__m256i ymm_V_dl16l_subG = _mm256_mullo_epi16(ymm_V0_16l, ymm_wKgv);
					__m256i ymm_V_dl16h_subG = _mm256_mullo_epi16(ymm_V0_16h, ymm_wKgv);

					__m256i ymm_V_dh16l_subG = _mm256_mullo_epi16(ymm_V1_16l, ymm_wKgv);
					__m256i ymm_V_dh16h_subG = _mm256_mullo_epi16(ymm_V1_16h, ymm_wKgv);


					__m256i ymm_U_dl16l_subG = _mm256_mullo_epi16(ymm_U0_16l, ymm_wKgu);
					__m256i ymm_U_dl16h_subG = _mm256_mullo_epi16(ymm_U0_16h, ymm_wKgu);

					__m256i ymm_U_dh16l_subG = _mm256_mullo_epi16(ymm_U1_16l, ymm_wKgu);
					__m256i ymm_U_dh16h_subG = _mm256_mullo_epi16(ymm_U1_16h, ymm_wKgu);

					ymm_V_dl16l_subG = _mm256_add_epi16(ymm_V_dl16l_subG, ymm_w_rounder);
					ymm_V_dl16h_subG = _mm256_add_epi16(ymm_V_dl16h_subG, ymm_w_rounder);

					ymm_V_dh16l_subG = _mm256_add_epi16(ymm_V_dh16l_subG, ymm_w_rounder);
					ymm_V_dh16h_subG = _mm256_add_epi16(ymm_V_dh16h_subG, ymm_w_rounder);

					ymm_U_dl16l_subG = _mm256_add_epi16(ymm_U_dl16l_subG, ymm_w_rounder);
					ymm_U_dl16h_subG = _mm256_add_epi16(ymm_U_dl16h_subG, ymm_w_rounder);

					ymm_U_dh16l_subG = _mm256_add_epi16(ymm_U_dh16l_subG, ymm_w_rounder);
					ymm_U_dh16h_subG = _mm256_add_epi16(ymm_U_dh16h_subG, ymm_w_rounder);


					ymm_V_dl16l_subG = _mm256_srai_epi16(ymm_V_dl16l_subG, 6);
					ymm_V_dl16h_subG = _mm256_srai_epi16(ymm_V_dl16h_subG, 6);

					ymm_V_dh16l_subG = _mm256_srai_epi16(ymm_V_dh16l_subG, 6);
					ymm_V_dh16h_subG = _mm256_srai_epi16(ymm_V_dh16h_subG, 6);


					ymm_U_dl16l_subG = _mm256_srai_epi16(ymm_U_dl16l_subG, 6);
					ymm_U_dl16h_subG = _mm256_srai_epi16(ymm_U_dl16h_subG, 6);

					ymm_U_dh16l_subG = _mm256_srai_epi16(ymm_U_dh16l_subG, 6);
					ymm_U_dh16h_subG = _mm256_srai_epi16(ymm_U_dh16h_subG, 6);


					__m256i ymm_G_0_16l = _mm256_sub_epi16(ymm_Y0_16l, ymm_V_dl16l_subG);
					__m256i ymm_G_0_16h = _mm256_sub_epi16(ymm_Y0_16h, ymm_V_dl16h_subG);

					ymm_G_0_16l = _mm256_sub_epi16(ymm_G_0_16l, ymm_U_dl16l_subG);
					ymm_G_0_16h = _mm256_sub_epi16(ymm_G_0_16h, ymm_U_dl16h_subG);

					__m256i ymm_G_1_16l = _mm256_sub_epi16(ymm_Y1_16l, ymm_V_dh16l_subG);
					__m256i ymm_G_1_16h = _mm256_sub_epi16(ymm_Y1_16h, ymm_V_dh16h_subG);

					ymm_G_1_16l = _mm256_sub_epi16(ymm_G_1_16l, ymm_U_dh16l_subG);
					ymm_G_1_16h = _mm256_sub_epi16(ymm_G_1_16h, ymm_U_dh16h_subG);

					// RGB post processing with gain and offset 
					__m256i ymm_RGBoffset = _mm256_set1_epi16(RGBo);
					__m256i ymm_RGBgain = _mm256_set1_epi16(RGBg);

					ymm_R_0_16l = _mm256_add_epi16(ymm_R_0_16l, ymm_RGBoffset);
					ymm_R_0_16h = _mm256_add_epi16(ymm_R_0_16h, ymm_RGBoffset);
					ymm_R_1_16l = _mm256_add_epi16(ymm_R_1_16l, ymm_RGBoffset);
					ymm_R_1_16h = _mm256_add_epi16(ymm_R_1_16h, ymm_RGBoffset);

					ymm_G_0_16l = _mm256_add_epi16(ymm_G_0_16l, ymm_RGBoffset);
					ymm_G_0_16h = _mm256_add_epi16(ymm_G_0_16h, ymm_RGBoffset);
					ymm_G_1_16l = _mm256_add_epi16(ymm_G_1_16l, ymm_RGBoffset);
					ymm_G_1_16h = _mm256_add_epi16(ymm_G_1_16h, ymm_RGBoffset);

					ymm_B_0_16l = _mm256_add_epi16(ymm_B_0_16l, ymm_RGBoffset);
					ymm_B_0_16h = _mm256_add_epi16(ymm_B_0_16h, ymm_RGBoffset);
					ymm_B_1_16l = _mm256_add_epi16(ymm_B_1_16l, ymm_RGBoffset);
					ymm_B_1_16h = _mm256_add_epi16(ymm_B_1_16h, ymm_RGBoffset);


					ymm_R_0_16l = _mm256_mullo_epi16(ymm_R_0_16l, ymm_RGBgain);
					ymm_R_0_16h = _mm256_mullo_epi16(ymm_R_0_16h, ymm_RGBgain);
					ymm_R_1_16l = _mm256_mullo_epi16(ymm_R_1_16l, ymm_RGBgain);
					ymm_R_1_16h = _mm256_mullo_epi16(ymm_R_1_16h, ymm_RGBgain);

					ymm_G_0_16l = _mm256_mullo_epi16(ymm_G_0_16l, ymm_RGBgain);
					ymm_G_0_16h = _mm256_mullo_epi16(ymm_G_0_16h, ymm_RGBgain);
					ymm_G_1_16l = _mm256_mullo_epi16(ymm_G_1_16l, ymm_RGBgain);
					ymm_G_1_16h = _mm256_mullo_epi16(ymm_G_1_16h, ymm_RGBgain);

					ymm_B_0_16l = _mm256_mullo_epi16(ymm_B_0_16l, ymm_RGBgain);
					ymm_B_0_16h = _mm256_mullo_epi16(ymm_B_0_16h, ymm_RGBgain);
					ymm_B_1_16l = _mm256_mullo_epi16(ymm_B_1_16l, ymm_RGBgain);
					ymm_B_1_16h = _mm256_mullo_epi16(ymm_B_1_16h, ymm_RGBgain);


					ymm_R_0_16l = _mm256_add_epi16(ymm_R_0_16l, ymm_w_rounder);
					ymm_R_0_16h = _mm256_add_epi16(ymm_R_0_16h, ymm_w_rounder);
					ymm_R_1_16l = _mm256_add_epi16(ymm_R_1_16l, ymm_w_rounder);
					ymm_R_1_16h = _mm256_add_epi16(ymm_R_1_16h, ymm_w_rounder);

					ymm_G_0_16l = _mm256_add_epi16(ymm_G_0_16l, ymm_w_rounder);
					ymm_G_0_16h = _mm256_add_epi16(ymm_G_0_16h, ymm_w_rounder);
					ymm_G_1_16l = _mm256_add_epi16(ymm_G_1_16l, ymm_w_rounder);
					ymm_G_1_16h = _mm256_add_epi16(ymm_G_1_16h, ymm_w_rounder);

					ymm_B_0_16l = _mm256_add_epi16(ymm_B_0_16l, ymm_w_rounder);
					ymm_B_0_16h = _mm256_add_epi16(ymm_B_0_16h, ymm_w_rounder);
					ymm_B_1_16l = _mm256_add_epi16(ymm_B_1_16l, ymm_w_rounder);
					ymm_B_1_16h = _mm256_add_epi16(ymm_B_1_16h, ymm_w_rounder);


					ymm_R_0_16l = _mm256_srai_epi16(ymm_R_0_16l, RGB_DIV_SHIFT);
					ymm_R_0_16h = _mm256_srai_epi16(ymm_R_0_16h, RGB_DIV_SHIFT);
					ymm_R_1_16l = _mm256_srai_epi16(ymm_R_1_16l, RGB_DIV_SHIFT);
					ymm_R_1_16h = _mm256_srai_epi16(ymm_R_1_16h, RGB_DIV_SHIFT);

					ymm_G_0_16l = _mm256_srai_epi16(ymm_G_0_16l, RGB_DIV_SHIFT);
					ymm_G_0_16h = _mm256_srai_epi16(ymm_G_0_16h, RGB_DIV_SHIFT);
					ymm_G_1_16l = _mm256_srai_epi16(ymm_G_1_16l, RGB_DIV_SHIFT);
					ymm_G_1_16h = _mm256_srai_epi16(ymm_G_1_16h, RGB_DIV_SHIFT);

					ymm_B_0_16l = _mm256_srai_epi16(ymm_B_0_16l, RGB_DIV_SHIFT);
					ymm_B_0_16h = _mm256_srai_epi16(ymm_B_0_16h, RGB_DIV_SHIFT);
					ymm_B_1_16l = _mm256_srai_epi16(ymm_B_1_16l, RGB_DIV_SHIFT);
					ymm_B_1_16h = _mm256_srai_epi16(ymm_B_1_16h, RGB_DIV_SHIFT);


					//pack 16bit to 8bit
					__m256i ymm_R_0_8 = _mm256_packus_epi16(ymm_R_0_16l, ymm_R_0_16h);
					__m256i ymm_R_1_8 = _mm256_packus_epi16(ymm_R_1_16l, ymm_R_1_16h);

					__m256i ymm_G_0_8 = _mm256_packus_epi16(ymm_G_0_16l, ymm_G_0_16h);
					__m256i ymm_G_1_8 = _mm256_packus_epi16(ymm_G_1_16l, ymm_G_1_16h);

					__m256i ymm_B_0_8 = _mm256_packus_epi16(ymm_B_0_16l, ymm_B_0_16h);
					__m256i ymm_B_1_8 = _mm256_packus_epi16(ymm_B_1_16l, ymm_B_1_16h);

/*					_mm256_store_si256((__m256i*)l_dstp_R, ymm_R_0_8);
					_mm256_store_si256((__m256i*)(l_dstp_R + 32), ymm_R_1_8);

					_mm256_store_si256((__m256i*)l_dstp_G, ymm_G_0_8);
					_mm256_store_si256((__m256i*)(l_dstp_G + 32), ymm_G_1_8);

					_mm256_store_si256((__m256i*)l_dstp_B, ymm_B_0_8);
					_mm256_store_si256((__m256i*)(l_dstp_B + 32), ymm_B_1_8);
					*/

					// convert planar RGB to BGRA32

					__m256i ymm_blend_G = _mm256_setr_epi8(0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0);
					__m256i ymm_blend_R = _mm256_setr_epi8(0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0);

					__m256i ymm_shuf_B_0_7 = _mm256_setr_epi8(0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255, 255);
					__m256i ymm_shuf_B_8_15 = _mm256_setr_epi8(8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255, 255);

					__m256i ymm_shuf_G_0_7 = _mm256_setr_epi8(255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255);
					__m256i ymm_shuf_G_8_15 = _mm256_setr_epi8(255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255);

					__m256i ymm_shuf_R_0_7 = _mm256_setr_epi8(255, 255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255);
					__m256i ymm_shuf_R_8_15 = _mm256_setr_epi8(255, 255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255);


					__m256i ymm_R_0_15 = _mm256_permute4x64_epi64(ymm_R_0_8, 0x44);
					__m256i ymm_G_0_15 = _mm256_permute4x64_epi64(ymm_G_0_8, 0x44);
					__m256i ymm_B_0_15 = _mm256_permute4x64_epi64(ymm_B_0_8, 0x44);

					// B
					__m256i ymm_BGRA_0_7 = _mm256_shuffle_epi8(ymm_B_0_15, ymm_shuf_B_0_7);
					__m256i ymm_BGRA_8_15 = _mm256_shuffle_epi8(ymm_B_0_15, ymm_shuf_B_8_15);

					// G
					__m256i ymm_G_0_7 = _mm256_shuffle_epi8(ymm_G_0_15, ymm_shuf_G_0_7);
					__m256i ymm_G_8_15 = _mm256_shuffle_epi8(ymm_G_0_15, ymm_shuf_G_8_15);

					ymm_BGRA_0_7 = _mm256_blendv_epi8(ymm_BGRA_0_7, ymm_G_0_7, ymm_blend_G);
					ymm_BGRA_8_15 = _mm256_blendv_epi8(ymm_BGRA_8_15, ymm_G_8_15, ymm_blend_G);

					//R
					__m256i ymm_R_0_7 = _mm256_shuffle_epi8(ymm_R_0_15, ymm_shuf_R_0_7);
					__m256i ymm_R_8_15 = _mm256_shuffle_epi8(ymm_R_0_15, ymm_shuf_R_8_15);

					ymm_BGRA_0_7 = _mm256_blendv_epi8(ymm_BGRA_0_7, ymm_R_0_7, ymm_blend_R);
					ymm_BGRA_8_15 = _mm256_blendv_epi8(ymm_BGRA_8_15, ymm_R_8_15, ymm_blend_R);

					////// BGR 16_23 and 24_31
					__m256i ymm_R_16_31 = _mm256_permute4x64_epi64(ymm_R_0_8, 0xEE);
					__m256i ymm_G_16_31 = _mm256_permute4x64_epi64(ymm_G_0_8, 0xEE);
					__m256i ymm_B_16_31 = _mm256_permute4x64_epi64(ymm_B_0_8, 0xEE);

					// B
					__m256i ymm_BGRA_16_23 = _mm256_shuffle_epi8(ymm_B_16_31, ymm_shuf_B_0_7);
					__m256i ymm_BGRA_24_31 = _mm256_shuffle_epi8(ymm_B_16_31, ymm_shuf_B_8_15);

					// G
					__m256i ymm_G_16_23 = _mm256_shuffle_epi8(ymm_G_16_31, ymm_shuf_G_0_7);
					__m256i ymm_G_24_31 = _mm256_shuffle_epi8(ymm_G_16_31, ymm_shuf_G_8_15);

					ymm_BGRA_16_23 = _mm256_blendv_epi8(ymm_BGRA_16_23, ymm_G_16_23, ymm_blend_G);
					ymm_BGRA_24_31 = _mm256_blendv_epi8(ymm_BGRA_24_31, ymm_G_24_31, ymm_blend_G);

					//R
					__m256i ymm_R_16_23 = _mm256_shuffle_epi8(ymm_R_16_31, ymm_shuf_R_0_7);
					__m256i ymm_R_24_31 = _mm256_shuffle_epi8(ymm_R_16_31, ymm_shuf_R_8_15);

					ymm_BGRA_16_23 = _mm256_blendv_epi8(ymm_BGRA_16_23, ymm_R_16_23, ymm_blend_R);
					ymm_BGRA_24_31 = _mm256_blendv_epi8(ymm_BGRA_24_31, ymm_R_24_31, ymm_blend_R);

					////// BGR 32_39 and 40_47
					__m256i ymm_R_32_47 = _mm256_permute4x64_epi64(ymm_R_1_8, 0x44);
					__m256i ymm_G_32_47 = _mm256_permute4x64_epi64(ymm_G_1_8, 0x44);
					__m256i ymm_B_32_47 = _mm256_permute4x64_epi64(ymm_B_1_8, 0x44);

					// B
					__m256i ymm_BGRA_32_39 = _mm256_shuffle_epi8(ymm_B_32_47, ymm_shuf_B_0_7);
					__m256i ymm_BGRA_40_47 = _mm256_shuffle_epi8(ymm_B_32_47, ymm_shuf_B_8_15);

					// G
					__m256i ymm_G_32_39 = _mm256_shuffle_epi8(ymm_G_32_47, ymm_shuf_G_0_7);
					__m256i ymm_G_40_47 = _mm256_shuffle_epi8(ymm_G_32_47, ymm_shuf_G_8_15);

					ymm_BGRA_32_39 = _mm256_blendv_epi8(ymm_BGRA_32_39, ymm_G_32_39, ymm_blend_G);
					ymm_BGRA_40_47 = _mm256_blendv_epi8(ymm_BGRA_40_47, ymm_G_40_47, ymm_blend_G);

					//R
					__m256i ymm_R_32_39 = _mm256_shuffle_epi8(ymm_R_32_47, ymm_shuf_R_0_7);
					__m256i ymm_R_40_47 = _mm256_shuffle_epi8(ymm_R_32_47, ymm_shuf_R_8_15);

					ymm_BGRA_32_39 = _mm256_blendv_epi8(ymm_BGRA_32_39, ymm_R_32_39, ymm_blend_R);
					ymm_BGRA_40_47 = _mm256_blendv_epi8(ymm_BGRA_40_47, ymm_R_40_47, ymm_blend_R);

					////// BGR 48_55 and 56_63
					__m256i ymm_R_48_63 = _mm256_permute4x64_epi64(ymm_R_1_8, 0xEE);
					__m256i ymm_G_48_63 = _mm256_permute4x64_epi64(ymm_G_1_8, 0xEE);
					__m256i ymm_B_48_63 = _mm256_permute4x64_epi64(ymm_B_1_8, 0xEE);

					// B
					__m256i ymm_BGRA_48_55 = _mm256_shuffle_epi8(ymm_B_48_63, ymm_shuf_B_0_7);
					__m256i ymm_BGRA_56_63 = _mm256_shuffle_epi8(ymm_B_48_63, ymm_shuf_B_8_15);

					// G
					__m256i ymm_G_48_55 = _mm256_shuffle_epi8(ymm_G_48_63, ymm_shuf_G_0_7);
					__m256i ymm_G_56_63 = _mm256_shuffle_epi8(ymm_G_48_63, ymm_shuf_G_8_15);

					ymm_BGRA_48_55 = _mm256_blendv_epi8(ymm_BGRA_48_55, ymm_G_48_55, ymm_blend_G);
					ymm_BGRA_56_63 = _mm256_blendv_epi8(ymm_BGRA_56_63, ymm_G_56_63, ymm_blend_G);

					//R
					__m256i ymm_R_48_55 = _mm256_shuffle_epi8(ymm_R_48_63, ymm_shuf_R_0_7);
					__m256i ymm_R_56_63 = _mm256_shuffle_epi8(ymm_R_48_63, ymm_shuf_R_8_15);

					ymm_BGRA_48_55 = _mm256_blendv_epi8(ymm_BGRA_48_55, ymm_R_48_55, ymm_blend_R);
					ymm_BGRA_56_63 = _mm256_blendv_epi8(ymm_BGRA_56_63, ymm_R_56_63, ymm_blend_R);

					if (bCacheStore)
					{
						_mm256_store_si256((__m256i*)l_dstp_BGRA, ymm_BGRA_0_7);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 32), ymm_BGRA_8_15);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 64), ymm_BGRA_16_23);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 96), ymm_BGRA_24_31);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 128), ymm_BGRA_32_39);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 160), ymm_BGRA_40_47);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 192), ymm_BGRA_48_55);
						_mm256_store_si256((__m256i*)(l_dstp_BGRA + 224), ymm_BGRA_56_63);
					}
					else
					{
						_mm256_stream_si256((__m256i*)l_dstp_BGRA, ymm_BGRA_0_7);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 32), ymm_BGRA_8_15);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 64), ymm_BGRA_16_23);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 96), ymm_BGRA_24_31);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 128), ymm_BGRA_32_39);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 160), ymm_BGRA_40_47);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 192), ymm_BGRA_48_55);
						_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 224), ymm_BGRA_56_63);
					}

					if (cs == 1)
					{
						l_srcp_Y += 64; // in bytes
						l_srcp_U += 32;
						l_srcp_V += 32;
					}
					else
					{
						l_srcp_Y += 128; // in bytes
						l_srcp_U += 64;
						l_srcp_V += 64;
					}

					l_dstp_R += 64;
					l_dstp_G += 64;
					l_dstp_B += 64;

					l_dstp_BGRA += 64 * 4;
				}

				// last cols
				int iUVadv = 0;
				for (int col = col64; col < row_size_Ysamples; ++col)
				{
					int iY;
					int iU;
					int iV;

					if (cs == 1)
					{
						iY = *l_srcp_Y;
						iU = *l_srcp_U - 128;
						iV = *l_srcp_V - 128;
					}
					else
					{
						iY = (int)*((unsigned short*)l_srcp_Y);
						iU = (int)*((unsigned short*)l_srcp_U);
						iV = (int)*((unsigned short*)l_srcp_V);

						if (bps == 10)
						{
							iY = iY >> 2;
							iU = iU >> 2;
							iV = iV >> 2;
						}

						if (bps == 12)
						{
							iY = iY >> 4;
							iU = iU >> 4;
							iV = iV >> 4;
						}

						if (bps == 14)
						{
							iY = iY >> 6;
							iU = iU >> 6;
							iV = iV >> 6;
						}

						if (bps == 16)
						{
							iY = iY >> 8;
							iU = iU >> 8;
							iV = iV >> 8;
						}

						iU = iU - 128;
						iV = iV - 128;

					}

					int iR = iY + (((iV * Kr) + RND_16) >> 6);
					int iB = iY + (((iU * Kb) + RND_16) >> 6);
					int iG = iY - (((iU * Kgu) + RND_16) >> 6) - (((iV * Kgv) + RND_16) >> 6);

					iR = (((iR + RGBo) * RGBg) + RND_16) >> RGB_DIV_SHIFT;
					iG = (((iG + RGBo) * RGBg) + RND_16) >> RGB_DIV_SHIFT;
					iB = (((iB + RGBo) * RGBg) + RND_16) >> RGB_DIV_SHIFT;

					if (iR < 0) iR = 0; if (iR > 255) iR = 255;
					if (iG < 0) iG = 0; if (iG > 255) iG = 255;
					if (iB < 0) iB = 0; if (iB > 255) iB = 255;

					int iBGRA = 0 | (unsigned char)iR << 16 | (unsigned char)iG << 8 | (unsigned char)iB;
					*(int*)(l_dstp_BGRA) = iBGRA;


					if (cs == 1)
					{
						l_srcp_Y += 1; // in bytes
						if (iUVadv % 2 != 0)
						{
							l_srcp_U += 1;
							l_srcp_V += 1;
						}
					}
					else
					{
						l_srcp_Y += 2; // in bytes
						if (iUVadv % 2 != 0)
						{
							l_srcp_U += 2;
							l_srcp_V += 2;
						}
					}

					iUVadv++;
					l_dstp_BGRA += 4;

				}
				
			}

			
	}

	if (!bCacheStore)
		_mm_sfence();
}


template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
void DecodeYUVtoRGB::DecodeYUV420imm32(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags)
{
	auto srcp_Y = src->GetReadPtr(PLANAR_Y);
	auto srcp_U = src->GetReadPtr(PLANAR_U);
	auto srcp_V = src->GetReadPtr(PLANAR_V);


	auto height = src->GetHeight(PLANAR_Y);
	auto row_size_Y = src->GetRowSize(PLANAR_Y);
	auto row_size_U = src->GetRowSize(PLANAR_U);
	auto row_size_V = src->GetRowSize(PLANAR_V);
	auto src_pitch_Y = src->GetPitch(PLANAR_Y);
	auto src_pitch_U = src->GetPitch(PLANAR_U);
	auto src_pitch_V = src->GetPitch(PLANAR_V);

	auto dstp_R = dst->GetWritePtr(4);
	auto dstp_G = dst->GetWritePtr(6);
	auto dstp_B = dst->GetWritePtr(2);

	auto dstp_BGRA = dst->GetWritePtr(2);
	auto dst_pitch_BGRA = dst->GetPitch();

	auto dst_pitch_R = dst->GetPitch(PLANAR_R);
	auto dst_pitch_G = dst->GetPitch(PLANAR_G);
	auto dst_pitch_B = dst->GetPitch(PLANAR_B);


#pragma omp parallel for num_threads(threads)
	for (uint64_t y = 0; y < height; y++)
	{
		unsigned char* l_dstp_R = dstp_R + y * dst_pitch_R;
		unsigned char* l_dstp_G = dstp_G + y * dst_pitch_G;
		unsigned char* l_dstp_B = dstp_B + y * dst_pitch_B;

		unsigned char* l_dstp_BGRA = dstp_BGRA + (height - y - 1) * dst_pitch_BGRA; // reverse scan for RGB interleaved

		unsigned char* l_srcp_Y = (unsigned char*)srcp_Y + y * src_pitch_Y;
		unsigned char* l_srcp_U = (unsigned char*)srcp_U + (y >> 1)* src_pitch_U;
		unsigned char* l_srcp_V = (unsigned char*)srcp_V + (y >> 1)* src_pitch_V;


		if (cpuFlags & CPUF_AVX2) // use AVX2
		{
			int col32;
			int row_size_Ysamples;
			if (cs == 1) // 8bit
			{
				col32 = row_size_Y - (row_size_Y % 32);
				row_size_Ysamples = row_size_Y;
			}
			else // 10 to 16 bit
			{
				col32 = (row_size_Y / 2) - ((row_size_Y / 2) % 32);
				row_size_Ysamples = row_size_Y / 2;
			}


			/*				int row_proc_size;
							if (cs == 1) // 8bit
								row_proc_size = row_size_Y;
							else // 10 to 16 bit
								row_proc_size = (row_size_Y / 2);
								*/
			__m256i ymm_dw_cbias = _mm256_set1_epi32(128);
			__m256i ymm_dw_rounder = _mm256_set1_epi32(RND_32);

			__m256i ymm_dwKr = _mm256_set1_epi32(Kr);
			__m256i ymm_dwKb = _mm256_set1_epi32(Kb);

			__m256i ymm_dwKgu = _mm256_set1_epi32(Kgu);
			__m256i ymm_dwKgv = _mm256_set1_epi32(Kgv);

			int UVpref = 0;


			// fill Y with 0 to 32 - debug
			
/*			for (int idx = 0; idx < 128; idx++)
			{
				l_srcp_Y[idx] = (unsigned int)idx;
			}
	*/		
			/*				unsigned short* lus_srcp_Y = (unsigned short*)l_srcp_Y;
							for (int idx = 0; idx < 64; idx++)
							{
								lus_srcp_Y[idx] = (unsigned short)(idx * 4);
							}
				*/
				
				// fill U with 0 to 16 - debug
		/*		for (int idx = 0; idx < 64; idx++)
				{
					l_srcp_U[idx] = (unsigned char)(idx + 128);
					l_srcp_V[idx] = (unsigned char)(idx + 128);
				}
			*/	
				// fill V with //0 to 32 - debug
/*				for (int idx = 0; idx < 16; idx++)
				{
					l_srcp_U[idx] = 128;
					l_srcp_V[idx] = 128;
				}
	*/			
				/*				unsigned short* lus_srcp_U = (unsigned short*)l_srcp_U;
								for (int idx = 0; idx < 32; idx++)
								{
									lus_srcp_U[idx] = (unsigned short)(idx * 4);
								}
				*/

			for (int col = 0; col < col32; col += 32)
			{
				__m256i ymm_Y0_32l;
				__m256i ymm_Y0_32h;

				__m256i ymm_Y1_32h;
				__m256i ymm_Y1_32l;

				__m256i ymm_U0_32l;
				__m256i ymm_V0_32l;

				__m256i ymm_U0_32h;
				__m256i ymm_V0_32h;

				__m256i ymm_U1_32l;
				__m256i ymm_V1_32l;

				__m256i ymm_U1_32h;
				__m256i ymm_V1_32h;


				if (cs == 1)
				{

					if (!bCacheLoad)
					{
						if (UVpref % 2 == 0)
						{
							_mm_prefetch((const CHAR*)(l_srcp_Y + 64), _MM_HINT_NTA);

							if (UVpref % 4 == 0)
							{
								_mm_prefetch((const CHAR*)(l_srcp_U + 64), _MM_HINT_NTA);
								_mm_prefetch((const CHAR*)(l_srcp_V + 64), _MM_HINT_NTA);
							}
						}

						UVpref++;
					}

					__m256i ymm0_Y0 = _mm256_load_si256((const __m256i*)l_srcp_Y); // should always load from 64-byte aligned start of row in AVS+ 3.7.3 (and later ?)
					
					__m256i ymm_Y_16l = _mm256_unpacklo_epi8(ymm0_Y0, _mm256_setzero_si256());
					__m256i ymm_Y_16h = _mm256_unpackhi_epi8(ymm0_Y0, _mm256_setzero_si256());

					ymm_Y0_32l = _mm256_unpacklo_epi16(ymm_Y_16l, _mm256_setzero_si256());
					ymm_Y0_32h = _mm256_unpackhi_epi16(ymm_Y_16l, _mm256_setzero_si256());

					ymm_Y1_32l = _mm256_unpacklo_epi16(ymm_Y_16h, _mm256_setzero_si256());
					ymm_Y1_32h = _mm256_unpackhi_epi16(ymm_Y_16h, _mm256_setzero_si256());

					__m128i xmm_U = _mm_load_si128((const __m128i*)(l_srcp_U));
					__m256i ymm_U = _mm256_permute4x64_epi64(_mm256_castsi128_si256(xmm_U), 0x50);
					
					__m128i xmm_V = _mm_load_si128((const __m128i*)(l_srcp_V));
					__m256i ymm_V = _mm256_permute4x64_epi64(_mm256_castsi128_si256(xmm_V), 0x50);

					__m256i ymm_U_8dl = _mm256_unpacklo_epi8(ymm_U, ymm_U);
					__m256i ymm_U_8dh = _mm256_unpackhi_epi8(ymm_U, ymm_U);

					__m256i ymm_U_16l = _mm256_unpacklo_epi8(ymm_U_8dl, _mm256_setzero_si256());
					__m256i ymm_U_16h = _mm256_unpackhi_epi8(ymm_U_8dl, _mm256_setzero_si256());

					__m256i ymm_V_8dl = _mm256_unpacklo_epi8(ymm_V, ymm_V);
					__m256i ymm_V_8dh = _mm256_unpackhi_epi8(ymm_V, ymm_V);

					__m256i ymm_V_16l = _mm256_unpacklo_epi8(ymm_V_8dl, _mm256_setzero_si256());
					__m256i ymm_V_16h = _mm256_unpackhi_epi8(ymm_V_8dl, _mm256_setzero_si256());


					ymm_U0_32l = _mm256_unpacklo_epi16(ymm_U_16l, _mm256_setzero_si256());
					ymm_V0_32l = _mm256_unpacklo_epi16(ymm_V_16l, _mm256_setzero_si256());

					ymm_U0_32h = _mm256_unpackhi_epi16(ymm_U_16l, _mm256_setzero_si256());
					ymm_V0_32h = _mm256_unpackhi_epi16(ymm_V_16l, _mm256_setzero_si256());

					ymm_U1_32l = _mm256_unpacklo_epi16(ymm_U_16h, _mm256_setzero_si256());
					ymm_V1_32l = _mm256_unpacklo_epi16(ymm_V_16h, _mm256_setzero_si256());

					ymm_U1_32h = _mm256_unpackhi_epi16(ymm_U_16h, _mm256_setzero_si256());
					ymm_V1_32h = _mm256_unpackhi_epi16(ymm_V_16h, _mm256_setzero_si256());


				}
/*				if (cs == 2)
				{
					if (!bCacheLoad)
					{
						_mm_prefetch((const CHAR*)(l_srcp_Y + 128), _MM_HINT_NTA);
						_mm_prefetch((const CHAR*)(l_srcp_Y + 128 + 64), _MM_HINT_NTA);

						_mm_prefetch((const CHAR*)(l_srcp_U + 64), _MM_HINT_NTA);
						_mm_prefetch((const CHAR*)(l_srcp_V + 64), _MM_HINT_NTA);
					}

					__m256i ymm_Y0_l = _mm256_load_si256((const __m256i*)l_srcp_Y); // should always load from 64-bit aligned start of row in AVS+ 3.7.3 (and later ?)
					__m256i ymm_Y0_h = _mm256_load_si256((const __m256i*)(l_srcp_Y + 32));
					__m256i ymm_Y1_l = _mm256_load_si256((const __m256i*)(l_srcp_Y + 64));
					__m256i ymm_Y1_h = _mm256_load_si256((const __m256i*)(l_srcp_Y + 96));

					__m256i ymm_U_l = _mm256_load_si256((const __m256i*)(l_srcp_U));
					__m256i ymm_U_h = _mm256_load_si256((const __m256i*)(l_srcp_U + 32));

					__m256i ymm_V_l = _mm256_load_si256((const __m256i*)(l_srcp_V));
					__m256i ymm_V_h = _mm256_load_si256((const __m256i*)(l_srcp_V + 32));

					ymm_Y0_16l = _mm256_permute2x128_si256(ymm_Y0_l, ymm_Y0_h, 0x20);
					ymm_Y0_16h = _mm256_permute2x128_si256(ymm_Y0_l, ymm_Y0_h, 0x31);

					ymm_Y1_16l = _mm256_permute2x128_si256(ymm_Y1_l, ymm_Y1_h, 0x20);
					ymm_Y1_16h = _mm256_permute2x128_si256(ymm_Y1_l, ymm_Y1_h, 0x31);

					if (bps == 10)
					{
						ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 2);
						ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 2);

						ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 2);
						ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 2);

						ymm_U_l = _mm256_srli_epi16(ymm_U_l, 2);
						ymm_U_h = _mm256_srli_epi16(ymm_U_h, 2);

						ymm_V_l = _mm256_srli_epi16(ymm_V_l, 2);
						ymm_V_h = _mm256_srli_epi16(ymm_V_h, 2);
					}

					if (bps == 12)
					{
						ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 4);
						ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 4);

						ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 4);
						ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 4);

						ymm_U_l = _mm256_srli_epi16(ymm_U_l, 4);
						ymm_U_h = _mm256_srli_epi16(ymm_U_h, 4);

						ymm_V_l = _mm256_srli_epi16(ymm_V_l, 4);
						ymm_V_h = _mm256_srli_epi16(ymm_V_h, 4);
					}

					if (bps == 14)
					{
						ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 6);
						ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 6);

						ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 6);
						ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 6);

						ymm_U_l = _mm256_srli_epi16(ymm_U_l, 6);
						ymm_U_h = _mm256_srli_epi16(ymm_U_h, 6);

						ymm_V_l = _mm256_srli_epi16(ymm_V_l, 6);
						ymm_V_h = _mm256_srli_epi16(ymm_V_h, 6);
					}

					if (bps == 16)
					{
						ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 8);
						ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 8);

						ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 8);
						ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 8);

						ymm_U_l = _mm256_srli_epi16(ymm_U_l, 8);
						ymm_U_h = _mm256_srli_epi16(ymm_U_h, 8);

						ymm_V_l = _mm256_srli_epi16(ymm_V_l, 8);
						ymm_V_h = _mm256_srli_epi16(ymm_V_h, 8);
					}


					ymm_U0_16l = _mm256_unpacklo_epi16(ymm_U_l, ymm_U_l);
					ymm_V0_16l = _mm256_unpacklo_epi16(ymm_V_l, ymm_V_l);

					ymm_U1_16l = _mm256_unpacklo_epi16(ymm_U_h, ymm_U_h);
					ymm_V1_16l = _mm256_unpacklo_epi16(ymm_V_h, ymm_V_h);

					ymm_U0_16h = _mm256_unpackhi_epi16(ymm_U_l, ymm_U_l);
					ymm_V0_16h = _mm256_unpackhi_epi16(ymm_V_l, ymm_V_l);

					ymm_U1_16h = _mm256_unpackhi_epi16(ymm_U_h, ymm_U_h);
					ymm_V1_16h = _mm256_unpackhi_epi16(ymm_V_h, ymm_V_h);


				}
				*/
				ymm_U0_32l = _mm256_sub_epi32(ymm_U0_32l, ymm_dw_cbias); // superscalarity of 3 at IceLake/SkyLake
				ymm_V0_32l = _mm256_sub_epi32(ymm_V0_32l, ymm_dw_cbias);

				ymm_U0_32h = _mm256_sub_epi32(ymm_U0_32h, ymm_dw_cbias);
				ymm_V0_32h = _mm256_sub_epi32(ymm_V0_32h, ymm_dw_cbias);

				ymm_U1_32l = _mm256_sub_epi32(ymm_U1_32l, ymm_dw_cbias); // superscalarity of 3 at IceLake/SkyLake
				ymm_V1_32l = _mm256_sub_epi32(ymm_V1_32l, ymm_dw_cbias);

				ymm_U1_32h = _mm256_sub_epi32(ymm_U1_32h, ymm_dw_cbias);
				ymm_V1_32h = _mm256_sub_epi32(ymm_V1_32h, ymm_dw_cbias);


				__m256i ymm_V_dl32l_addR = _mm256_mullo_epi32(ymm_V0_32l, ymm_dwKr);
				__m256i ymm_V_dl32h_addR = _mm256_mullo_epi32(ymm_V0_32h, ymm_dwKr);

				__m256i ymm_V_dh32l_addR = _mm256_mullo_epi32(ymm_V1_32l, ymm_dwKr);
				__m256i ymm_V_dh32h_addR = _mm256_mullo_epi32(ymm_V1_32h, ymm_dwKr);


				__m256i ymm_U_dl32l_addB = _mm256_mullo_epi32(ymm_U0_32l, ymm_dwKb);
				__m256i ymm_U_dl32h_addB = _mm256_mullo_epi32(ymm_U0_32h, ymm_dwKb);

				__m256i ymm_U_dh32l_addB = _mm256_mullo_epi32(ymm_U1_32l, ymm_dwKb);
				__m256i ymm_U_dh32h_addB = _mm256_mullo_epi32(ymm_U1_32h, ymm_dwKb);


				ymm_V_dl32l_addR = _mm256_add_epi32(ymm_V_dl32l_addR, ymm_dw_rounder);
				ymm_V_dl32h_addR = _mm256_add_epi32(ymm_V_dl32h_addR, ymm_dw_rounder);

				ymm_V_dh32l_addR = _mm256_add_epi32(ymm_V_dh32l_addR, ymm_dw_rounder);
				ymm_V_dh32h_addR = _mm256_add_epi32(ymm_V_dh32h_addR, ymm_dw_rounder);

				ymm_U_dl32l_addB = _mm256_add_epi32(ymm_U_dl32l_addB, ymm_dw_rounder);
				ymm_U_dl32h_addB = _mm256_add_epi32(ymm_U_dl32h_addB, ymm_dw_rounder);

				ymm_U_dh32l_addB = _mm256_add_epi32(ymm_U_dh32l_addB, ymm_dw_rounder);
				ymm_U_dh32h_addB = _mm256_add_epi32(ymm_U_dh32h_addB, ymm_dw_rounder);


				ymm_V_dl32l_addR = _mm256_srai_epi32(ymm_V_dl32l_addR, 13);
				ymm_V_dl32h_addR = _mm256_srai_epi32(ymm_V_dl32h_addR, 13);

				ymm_V_dh32l_addR = _mm256_srai_epi32(ymm_V_dh32l_addR, 13);
				ymm_V_dh32h_addR = _mm256_srai_epi32(ymm_V_dh32h_addR, 13);

				ymm_U_dl32l_addB = _mm256_srai_epi32(ymm_U_dl32l_addB, 13);
				ymm_U_dl32h_addB = _mm256_srai_epi32(ymm_U_dl32h_addB, 13);

				ymm_U_dh32l_addB = _mm256_srai_epi32(ymm_U_dh32l_addB, 13);
				ymm_U_dh32h_addB = _mm256_srai_epi32(ymm_U_dh32h_addB, 13);


				// R
				__m256i ymm_R_0_32l = _mm256_add_epi32(ymm_Y0_32l, ymm_V_dl32l_addR);
				__m256i ymm_R_0_32h = _mm256_add_epi32(ymm_Y0_32h, ymm_V_dl32h_addR);

				__m256i ymm_R_1_32l = _mm256_add_epi32(ymm_Y1_32l, ymm_V_dh32l_addR);
				__m256i ymm_R_1_32h = _mm256_add_epi32(ymm_Y1_32h, ymm_V_dh32h_addR);

				// B
				__m256i ymm_B_0_32l = _mm256_add_epi32(ymm_Y0_32l, ymm_U_dl32l_addB);
				__m256i ymm_B_0_32h = _mm256_add_epi32(ymm_Y0_32h, ymm_U_dl32h_addB);

				__m256i ymm_B_1_32l = _mm256_add_epi32(ymm_Y1_32l, ymm_U_dh32l_addB);
				__m256i ymm_B_1_32h = _mm256_add_epi32(ymm_Y1_32h, ymm_U_dh32h_addB);

				//G
				__m256i ymm_V_dl32l_subG = _mm256_mullo_epi32(ymm_V0_32l, ymm_dwKgv);
				__m256i ymm_V_dl32h_subG = _mm256_mullo_epi32(ymm_V0_32h, ymm_dwKgv);

				__m256i ymm_V_dh32l_subG = _mm256_mullo_epi32(ymm_V1_32l, ymm_dwKgv);
				__m256i ymm_V_dh32h_subG = _mm256_mullo_epi32(ymm_V1_32h, ymm_dwKgv);


				__m256i ymm_U_dl32l_subG = _mm256_mullo_epi32(ymm_U0_32l, ymm_dwKgu);
				__m256i ymm_U_dl32h_subG = _mm256_mullo_epi32(ymm_U0_32h, ymm_dwKgu);

				__m256i ymm_U_dh32l_subG = _mm256_mullo_epi32(ymm_U1_32l, ymm_dwKgu);
				__m256i ymm_U_dh32h_subG = _mm256_mullo_epi32(ymm_U1_32h, ymm_dwKgu);


				ymm_V_dl32l_subG = _mm256_add_epi32(ymm_V_dl32l_subG, ymm_dw_rounder);
				ymm_V_dl32h_subG = _mm256_add_epi32(ymm_V_dl32h_subG, ymm_dw_rounder);

				ymm_V_dh32l_subG = _mm256_add_epi32(ymm_V_dh32l_subG, ymm_dw_rounder);
				ymm_V_dh32h_subG = _mm256_add_epi32(ymm_V_dh32h_subG, ymm_dw_rounder);

				ymm_U_dl32l_subG = _mm256_add_epi32(ymm_U_dl32l_subG, ymm_dw_rounder);
				ymm_U_dl32h_subG = _mm256_add_epi32(ymm_U_dl32h_subG, ymm_dw_rounder);

				ymm_U_dh32l_subG = _mm256_add_epi32(ymm_U_dh32l_subG, ymm_dw_rounder);
				ymm_U_dh32h_subG = _mm256_add_epi32(ymm_U_dh32h_subG, ymm_dw_rounder);


				ymm_V_dl32l_subG = _mm256_srai_epi32(ymm_V_dl32l_subG, 13);
				ymm_V_dl32h_subG = _mm256_srai_epi32(ymm_V_dl32h_subG, 13);

				ymm_V_dh32l_subG = _mm256_srai_epi32(ymm_V_dh32l_subG, 13);
				ymm_V_dh32h_subG = _mm256_srai_epi32(ymm_V_dh32h_subG, 13);

				ymm_U_dl32l_subG = _mm256_srai_epi32(ymm_U_dl32l_subG, 13);
				ymm_U_dl32h_subG = _mm256_srai_epi32(ymm_U_dl32h_subG, 13);

				ymm_U_dh32l_subG = _mm256_srai_epi32(ymm_U_dh32l_subG, 13);
				ymm_U_dh32h_subG = _mm256_srai_epi32(ymm_U_dh32h_subG, 13);


				__m256i ymm_G_0_32l = _mm256_sub_epi32(ymm_Y0_32l, ymm_V_dl32l_subG);
				__m256i ymm_G_0_32h = _mm256_sub_epi32(ymm_Y0_32h, ymm_V_dl32h_subG);

				ymm_G_0_32l = _mm256_sub_epi32(ymm_G_0_32l, ymm_U_dl32l_subG);
				ymm_G_0_32h = _mm256_sub_epi32(ymm_G_0_32h, ymm_U_dl32h_subG);

				__m256i ymm_G_1_32l = _mm256_sub_epi32(ymm_Y1_32l, ymm_V_dh32l_subG);
				__m256i ymm_G_1_32h = _mm256_sub_epi32(ymm_Y1_32h, ymm_V_dh32h_subG);

				ymm_G_1_32l = _mm256_sub_epi32(ymm_G_1_32l, ymm_U_dh32l_subG);
				ymm_G_1_32h = _mm256_sub_epi32(ymm_G_1_32h, ymm_U_dh32h_subG);

				// RGB post processing with gain and offset 
				__m256i ymm_RGBoffset = _mm256_set1_epi32(RGBo);
				__m256i ymm_RGBgain = _mm256_set1_epi32(RGBg);

				ymm_R_0_32l = _mm256_add_epi32(ymm_R_0_32l, ymm_RGBoffset);
				ymm_R_0_32h = _mm256_add_epi32(ymm_R_0_32h, ymm_RGBoffset);
				ymm_R_1_32l = _mm256_add_epi32(ymm_R_1_32l, ymm_RGBoffset);
				ymm_R_1_32h = _mm256_add_epi32(ymm_R_1_32h, ymm_RGBoffset);

				ymm_G_0_32l = _mm256_add_epi32(ymm_G_0_32l, ymm_RGBoffset);
				ymm_G_0_32h = _mm256_add_epi32(ymm_G_0_32h, ymm_RGBoffset);
				ymm_G_1_32l = _mm256_add_epi32(ymm_G_1_32l, ymm_RGBoffset);
				ymm_G_1_32h = _mm256_add_epi32(ymm_G_1_32h, ymm_RGBoffset);

				ymm_B_0_32l = _mm256_add_epi32(ymm_B_0_32l, ymm_RGBoffset);
				ymm_B_0_32h = _mm256_add_epi32(ymm_B_0_32h, ymm_RGBoffset);
				ymm_B_1_32l = _mm256_add_epi32(ymm_B_1_32l, ymm_RGBoffset);
				ymm_B_1_32h = _mm256_add_epi32(ymm_B_1_32h, ymm_RGBoffset);


				ymm_R_0_32l = _mm256_mullo_epi32(ymm_R_0_32l, ymm_RGBgain);
				ymm_R_0_32h = _mm256_mullo_epi32(ymm_R_0_32h, ymm_RGBgain);
				ymm_R_1_32l = _mm256_mullo_epi32(ymm_R_1_32l, ymm_RGBgain);
				ymm_R_1_32h = _mm256_mullo_epi32(ymm_R_1_32h, ymm_RGBgain);

				ymm_G_0_32l = _mm256_mullo_epi32(ymm_G_0_32l, ymm_RGBgain);
				ymm_G_0_32h = _mm256_mullo_epi32(ymm_G_0_32h, ymm_RGBgain);
				ymm_G_1_32l = _mm256_mullo_epi32(ymm_G_1_32l, ymm_RGBgain);
				ymm_G_1_32h = _mm256_mullo_epi32(ymm_G_1_32h, ymm_RGBgain);

				ymm_B_0_32l = _mm256_mullo_epi32(ymm_B_0_32l, ymm_RGBgain);
				ymm_B_0_32h = _mm256_mullo_epi32(ymm_B_0_32h, ymm_RGBgain);
				ymm_B_1_32l = _mm256_mullo_epi32(ymm_B_1_32l, ymm_RGBgain);
				ymm_B_1_32h = _mm256_mullo_epi32(ymm_B_1_32h, ymm_RGBgain);


				ymm_R_0_32l = _mm256_add_epi32(ymm_R_0_32l, ymm_dw_rounder);
				ymm_R_0_32h = _mm256_add_epi32(ymm_R_0_32h, ymm_dw_rounder);
				ymm_R_1_32l = _mm256_add_epi32(ymm_R_1_32l, ymm_dw_rounder);
				ymm_R_1_32h = _mm256_add_epi32(ymm_R_1_32h, ymm_dw_rounder);

				ymm_G_0_32l = _mm256_add_epi32(ymm_G_0_32l, ymm_dw_rounder);
				ymm_G_0_32h = _mm256_add_epi32(ymm_G_0_32h, ymm_dw_rounder);
				ymm_G_1_32l = _mm256_add_epi32(ymm_G_1_32l, ymm_dw_rounder);
				ymm_G_1_32h = _mm256_add_epi32(ymm_G_1_32h, ymm_dw_rounder);

				ymm_B_0_32l = _mm256_add_epi32(ymm_B_0_32l, ymm_dw_rounder);
				ymm_B_0_32h = _mm256_add_epi32(ymm_B_0_32h, ymm_dw_rounder);
				ymm_B_1_32l = _mm256_add_epi32(ymm_B_1_32l, ymm_dw_rounder);
				ymm_B_1_32h = _mm256_add_epi32(ymm_B_1_32h, ymm_dw_rounder);


				ymm_R_0_32l = _mm256_srai_epi32(ymm_R_0_32l, RGB_DIV_SHIFT_32);
				ymm_R_0_32h = _mm256_srai_epi32(ymm_R_0_32h, RGB_DIV_SHIFT_32);
				ymm_R_1_32l = _mm256_srai_epi32(ymm_R_1_32l, RGB_DIV_SHIFT_32);
				ymm_R_1_32h = _mm256_srai_epi32(ymm_R_1_32h, RGB_DIV_SHIFT_32);

				ymm_G_0_32l = _mm256_srai_epi32(ymm_G_0_32l, RGB_DIV_SHIFT_32);
				ymm_G_0_32h = _mm256_srai_epi32(ymm_G_0_32h, RGB_DIV_SHIFT_32);
				ymm_G_1_32l = _mm256_srai_epi32(ymm_G_1_32l, RGB_DIV_SHIFT_32);
				ymm_G_1_32h = _mm256_srai_epi32(ymm_G_1_32h, RGB_DIV_SHIFT_32);

				ymm_B_0_32l = _mm256_srai_epi32(ymm_B_0_32l, RGB_DIV_SHIFT_32);
				ymm_B_0_32h = _mm256_srai_epi32(ymm_B_0_32h, RGB_DIV_SHIFT_32);
				ymm_B_1_32l = _mm256_srai_epi32(ymm_B_1_32l, RGB_DIV_SHIFT_32);
				ymm_B_1_32h = _mm256_srai_epi32(ymm_B_1_32h, RGB_DIV_SHIFT_32);


				//pack 32bit to 8bit
				__m256i ymm_R_0_16 = _mm256_packus_epi32(ymm_R_0_32l, ymm_R_0_32h);
				__m256i ymm_R_1_16 = _mm256_packus_epi32(ymm_R_1_32l, ymm_R_1_32h);

				__m256i ymm_R_8 = _mm256_packus_epi16(ymm_R_0_16, ymm_R_1_16);

				__m256i ymm_G_0_16 = _mm256_packus_epi32(ymm_G_0_32l, ymm_G_0_32h);
				__m256i ymm_G_1_16 = _mm256_packus_epi32(ymm_G_1_32l, ymm_G_1_32h);

				__m256i ymm_G_8 = _mm256_packus_epi16(ymm_G_0_16, ymm_G_1_16);

				__m256i ymm_B_0_16 = _mm256_packus_epi32(ymm_B_0_32l, ymm_B_0_32h);
				__m256i ymm_B_1_16 = _mm256_packus_epi32(ymm_B_1_32l, ymm_B_1_32h);

				__m256i ymm_B_8 = _mm256_packus_epi16(ymm_B_0_16, ymm_B_1_16);

				/*
				if (bCacheStore)
				{
					_mm256_store_si256((__m256i*)l_dstp_R, ymm_R_8);
					_mm256_store_si256((__m256i*)l_dstp_G, ymm_G_8);
					_mm256_store_si256((__m256i*)l_dstp_B, ymm_B_8);
				}
				else
				{
					_mm256_stream_si256((__m256i*)l_dstp_R, ymm_R_8);
					_mm256_stream_si256((__m256i*)l_dstp_G, ymm_G_8);
					_mm256_stream_si256((__m256i*)l_dstp_B, ymm_B_8);

				}
				*/

				// convert planar RGB to BGRA32
				
				__m256i ymm_blend_G = _mm256_setr_epi8(0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0);
				__m256i ymm_blend_R = _mm256_setr_epi8(0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0);

				__m256i ymm_shuf_B_0_7 = _mm256_setr_epi8(0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255, 255);
				__m256i ymm_shuf_B_8_15 = _mm256_setr_epi8(8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255, 255);

				__m256i ymm_shuf_G_0_7 = _mm256_setr_epi8(255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255);
				__m256i ymm_shuf_G_8_15 = _mm256_setr_epi8(255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255);

				__m256i ymm_shuf_R_0_7 = _mm256_setr_epi8(255, 255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255);
				__m256i ymm_shuf_R_8_15 = _mm256_setr_epi8(255, 255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255);


				__m256i ymm_R_0_15 = _mm256_permute4x64_epi64(ymm_R_8, 0x44);
				__m256i ymm_G_0_15 = _mm256_permute4x64_epi64(ymm_G_8, 0x44);
				__m256i ymm_B_0_15 = _mm256_permute4x64_epi64(ymm_B_8, 0x44);

				// B
				__m256i ymm_BGRA_0_7 = _mm256_shuffle_epi8(ymm_B_0_15, ymm_shuf_B_0_7);
				__m256i ymm_BGRA_8_15 = _mm256_shuffle_epi8(ymm_B_0_15, ymm_shuf_B_8_15);

				// G
				__m256i ymm_G_0_7 = _mm256_shuffle_epi8(ymm_G_0_15, ymm_shuf_G_0_7);
				__m256i ymm_G_8_15 = _mm256_shuffle_epi8(ymm_G_0_15, ymm_shuf_G_8_15);

				ymm_BGRA_0_7 = _mm256_blendv_epi8(ymm_BGRA_0_7, ymm_G_0_7, ymm_blend_G);
				ymm_BGRA_8_15 = _mm256_blendv_epi8(ymm_BGRA_8_15, ymm_G_8_15, ymm_blend_G);

				//R
				__m256i ymm_R_0_7 = _mm256_shuffle_epi8(ymm_R_0_15, ymm_shuf_R_0_7);
				__m256i ymm_R_8_15 = _mm256_shuffle_epi8(ymm_R_0_15, ymm_shuf_R_8_15);

				ymm_BGRA_0_7 = _mm256_blendv_epi8(ymm_BGRA_0_7, ymm_R_0_7, ymm_blend_R);
				ymm_BGRA_8_15 = _mm256_blendv_epi8(ymm_BGRA_8_15, ymm_R_8_15, ymm_blend_R);

				////// BGR 16_23 and 24_31
				__m256i ymm_R_16_31 = _mm256_permute4x64_epi64(ymm_R_8, 0xEE);
				__m256i ymm_G_16_31 = _mm256_permute4x64_epi64(ymm_G_8, 0xEE);
				__m256i ymm_B_16_31 = _mm256_permute4x64_epi64(ymm_B_8, 0xEE);

				// B
				__m256i ymm_BGRA_16_23 = _mm256_shuffle_epi8(ymm_B_16_31, ymm_shuf_B_0_7);
				__m256i ymm_BGRA_24_31 = _mm256_shuffle_epi8(ymm_B_16_31, ymm_shuf_B_8_15);

				// G
				__m256i ymm_G_16_23 = _mm256_shuffle_epi8(ymm_G_16_31, ymm_shuf_G_0_7);
				__m256i ymm_G_24_31 = _mm256_shuffle_epi8(ymm_G_16_31, ymm_shuf_G_8_15);

				ymm_BGRA_16_23 = _mm256_blendv_epi8(ymm_BGRA_16_23, ymm_G_16_23, ymm_blend_G);
				ymm_BGRA_24_31 = _mm256_blendv_epi8(ymm_BGRA_24_31, ymm_G_24_31, ymm_blend_G);

				//R
				__m256i ymm_R_16_23 = _mm256_shuffle_epi8(ymm_R_16_31, ymm_shuf_R_0_7);
				__m256i ymm_R_24_31 = _mm256_shuffle_epi8(ymm_R_16_31, ymm_shuf_R_8_15);

				ymm_BGRA_16_23 = _mm256_blendv_epi8(ymm_BGRA_16_23, ymm_R_16_23, ymm_blend_R);
				ymm_BGRA_24_31 = _mm256_blendv_epi8(ymm_BGRA_24_31, ymm_R_24_31, ymm_blend_R);

				if (bCacheStore)
				{
					_mm256_store_si256((__m256i*)l_dstp_BGRA, ymm_BGRA_0_7);
					_mm256_store_si256((__m256i*)(l_dstp_BGRA + 32), ymm_BGRA_8_15);
					_mm256_store_si256((__m256i*)(l_dstp_BGRA + 64), ymm_BGRA_16_23);
					_mm256_store_si256((__m256i*)(l_dstp_BGRA + 96), ymm_BGRA_24_31);
				}
				else
				{
					_mm256_stream_si256((__m256i*)l_dstp_BGRA, ymm_BGRA_0_7);
					_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 32), ymm_BGRA_8_15);
					_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 64), ymm_BGRA_16_23);
					_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 96), ymm_BGRA_24_31);
				}
				
				if (cs == 1)
				{
					l_srcp_Y += 32; // in bytes
					l_srcp_U += 16;
					l_srcp_V += 16;
				}
				else
				{
					l_srcp_Y += 64; // in bytes
					l_srcp_U += 32;
					l_srcp_V += 32;
				}

				l_dstp_R += 32;
				l_dstp_G += 32;
				l_dstp_B += 32;

				l_dstp_BGRA += 32 * 4;
			}

			// last cols
			int iUVadv = 0;
			for (int col = col32; col < row_size_Ysamples; ++col)
			{
				int iY;
				int iU;
				int iV;

				if (cs == 1)
				{
					iY = *l_srcp_Y;
					iU = *l_srcp_U - 128;
					iV = *l_srcp_V - 128;
				}
				else
				{
					iY = (int) * ((unsigned short*)l_srcp_Y);
					iU = (int) * ((unsigned short*)l_srcp_U);
					iV = (int) * ((unsigned short*)l_srcp_V);

					if (bps == 10)
					{
						iY = iY >> 2;
						iU = iU >> 2;
						iV = iV >> 2;
					}

					if (bps == 12)
					{
						iY = iY >> 4;
						iU = iU >> 4;
						iV = iV >> 4;
					}

					if (bps == 14)
					{
						iY = iY >> 6;
						iU = iU >> 6;
						iV = iV >> 6;
					}

					if (bps == 16)
					{
						iY = iY >> 8;
						iU = iU >> 8;
						iV = iV >> 8;
					}

					iU = iU - 128;
					iV = iV - 128;

				}

				int iR = iY + (((iV * Kr) + RND_32) >> 13);
				int iB = iY + (((iU * Kb) + RND_32) >> 13);
				int iG = iY - (((iU * Kgu) + RND_32) >> 13) - (((iV * Kgv) + RND_32) >> 13);

				iR = (((iR + RGBo) * RGBg) + RND_32) >> RGB_DIV_SHIFT_32;
				iG = (((iG + RGBo) * RGBg) + RND_32) >> RGB_DIV_SHIFT_32;
				iB = (((iB + RGBo) * RGBg) + RND_32) >> RGB_DIV_SHIFT_32;

				if (iR < 0) iR = 0; if (iR > 255) iR = 255;
				if (iG < 0) iG = 0; if (iG > 255) iG = 255;
				if (iB < 0) iB = 0; if (iB > 255) iB = 255;

				int iBGRA = 0 | (unsigned char)iR << 16 | (unsigned char)iG << 8 | (unsigned char)iB;
				*(int*)(l_dstp_BGRA) = iBGRA;


				if (cs == 1)
				{
					l_srcp_Y += 1; // in bytes
					if (iUVadv % 2 != 0)
					{
						l_srcp_U += 1;
						l_srcp_V += 1;
					}
				}
				else
				{
					l_srcp_Y += 2; // in bytes
					if (iUVadv % 2 != 0)
					{
						l_srcp_U += 2;
						l_srcp_V += 2;
					}
				}

				iUVadv++;
				l_dstp_BGRA += 4;

			}

		}


	}

	if (!bCacheStore)
		_mm_sfence();
}

template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
void DecodeYUVtoRGB::DecodeYUV420imm32_mat(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags)
{
	auto srcp_Y = src->GetReadPtr(PLANAR_Y);
	auto srcp_U = src->GetReadPtr(PLANAR_U);
	auto srcp_V = src->GetReadPtr(PLANAR_V);


	auto height = src->GetHeight(PLANAR_Y);
	auto row_size_Y = src->GetRowSize(PLANAR_Y);
	auto row_size_U = src->GetRowSize(PLANAR_U);
	auto row_size_V = src->GetRowSize(PLANAR_V);
	auto src_pitch_Y = src->GetPitch(PLANAR_Y);
	auto src_pitch_U = src->GetPitch(PLANAR_U);
	auto src_pitch_V = src->GetPitch(PLANAR_V);

	auto dstp_R = dst->GetWritePtr(4);
	auto dstp_G = dst->GetWritePtr(6);
	auto dstp_B = dst->GetWritePtr(2);

	auto dstp_BGRA = dst->GetWritePtr(2);
	auto dst_pitch_BGRA = dst->GetPitch();

	auto dst_pitch_R = dst->GetPitch(PLANAR_R);
	auto dst_pitch_G = dst->GetPitch(PLANAR_G);
	auto dst_pitch_B = dst->GetPitch(PLANAR_B);


#pragma omp parallel for num_threads(threads)
	for (uint64_t y = 0; y < height; y++)
	{
		unsigned char* l_dstp_R = dstp_R + y * dst_pitch_R;
		unsigned char* l_dstp_G = dstp_G + y * dst_pitch_G;
		unsigned char* l_dstp_B = dstp_B + y * dst_pitch_B;

		unsigned char* l_dstp_BGRA = dstp_BGRA + (height - y - 1) * dst_pitch_BGRA; // reverse scan for RGB interleaved

		unsigned char* l_srcp_Y = (unsigned char*)srcp_Y + y * src_pitch_Y;
		unsigned char* l_srcp_U = (unsigned char*)srcp_U + (y >> 1)* src_pitch_U;
		unsigned char* l_srcp_V = (unsigned char*)srcp_V + (y >> 1)* src_pitch_V;


		if (cpuFlags & CPUF_AVX2) // use AVX2
		{
			int col32;
			int row_size_Ysamples;
			if (cs == 1) // 8bit
			{
				col32 = row_size_Y - (row_size_Y % 32);
				row_size_Ysamples = row_size_Y;
			}
			else // 10 to 16 bit
			{
				col32 = (row_size_Y / 2) - ((row_size_Y / 2) % 32);
				row_size_Ysamples = row_size_Y / 2;
			}


			/*				int row_proc_size;
							if (cs == 1) // 8bit
								row_proc_size = row_size_Y;
							else // 10 to 16 bit
								row_proc_size = (row_size_Y / 2);
								*/
			__m256i ymm_dw_cbias = _mm256_set1_epi32(128);

			__m256i ymm_dwKr = _mm256_set1_epi32(Kr);
			__m256i ymm_dwKb = _mm256_set1_epi32(Kb);

			__m256i ymm_dwKgu = _mm256_set1_epi32(Kgu);
			__m256i ymm_dwKgv = _mm256_set1_epi32(Kgv);

			int UVpref = 0;


			// fill Y with 0 to 32 - debug

/*			for (int idx = 0; idx < 128; idx++)
			{
				l_srcp_Y[idx] = (unsigned int)idx;
			}
	*/
	/*				unsigned short* lus_srcp_Y = (unsigned short*)l_srcp_Y;
					for (int idx = 0; idx < 64; idx++)
					{
						lus_srcp_Y[idx] = (unsigned short)(idx * 4);
					}
		*/

		// fill U with 0 to 16 - debug
/*		for (int idx = 0; idx < 64; idx++)
		{
			l_srcp_U[idx] = (unsigned char)(idx + 128);
			l_srcp_V[idx] = (unsigned char)(idx + 128);
		}
	*/
	// fill V with //0 to 32 - debug
/*				for (int idx = 0; idx < 16; idx++)
				{
					l_srcp_U[idx] = 128;
					l_srcp_V[idx] = 128;
				}
	*/
	/*				unsigned short* lus_srcp_U = (unsigned short*)l_srcp_U;
					for (int idx = 0; idx < 32; idx++)
					{
						lus_srcp_U[idx] = (unsigned short)(idx * 4);
					}
	*/

			for (int col = 0; col < col32; col += 32)
			{
				__m256i ymm_Y0_32l;
				__m256i ymm_Y0_32h;

				__m256i ymm_Y1_32h;
				__m256i ymm_Y1_32l;

				__m256i ymm_U0_32l;
				__m256i ymm_V0_32l;

				__m256i ymm_U0_32h;
				__m256i ymm_V0_32h;

				__m256i ymm_U1_32l;
				__m256i ymm_V1_32l;

				__m256i ymm_U1_32h;
				__m256i ymm_V1_32h;


				if (cs == 1)
				{

					if (!bCacheLoad)
					{
						if (UVpref % 2 == 0)
						{
							_mm_prefetch((const CHAR*)(l_srcp_Y + 64), _MM_HINT_NTA);

							if (UVpref % 4 == 0)
							{
								_mm_prefetch((const CHAR*)(l_srcp_U + 64), _MM_HINT_NTA);
								_mm_prefetch((const CHAR*)(l_srcp_V + 64), _MM_HINT_NTA);
							}
						}

						UVpref++;
					}

					__m256i ymm0_Y0 = _mm256_load_si256((const __m256i*)l_srcp_Y); // should always load from 64-byte aligned start of row in AVS+ 3.7.3 (and later ?)

					__m256i ymm_Y_16l = _mm256_unpacklo_epi8(ymm0_Y0, _mm256_setzero_si256());
					__m256i ymm_Y_16h = _mm256_unpackhi_epi8(ymm0_Y0, _mm256_setzero_si256());

					ymm_Y0_32l = _mm256_unpacklo_epi16(ymm_Y_16l, _mm256_setzero_si256());
					ymm_Y0_32h = _mm256_unpackhi_epi16(ymm_Y_16l, _mm256_setzero_si256());

					ymm_Y1_32l = _mm256_unpacklo_epi16(ymm_Y_16h, _mm256_setzero_si256());
					ymm_Y1_32h = _mm256_unpackhi_epi16(ymm_Y_16h, _mm256_setzero_si256());

					__m128i xmm_U = _mm_load_si128((const __m128i*)(l_srcp_U));
					__m256i ymm_U = _mm256_permute4x64_epi64(_mm256_castsi128_si256(xmm_U), 0x50);

					__m128i xmm_V = _mm_load_si128((const __m128i*)(l_srcp_V));
					__m256i ymm_V = _mm256_permute4x64_epi64(_mm256_castsi128_si256(xmm_V), 0x50);

					__m256i ymm_U_8dl = _mm256_unpacklo_epi8(ymm_U, ymm_U);
					__m256i ymm_U_8dh = _mm256_unpackhi_epi8(ymm_U, ymm_U);

					__m256i ymm_U_16l = _mm256_unpacklo_epi8(ymm_U_8dl, _mm256_setzero_si256());
					__m256i ymm_U_16h = _mm256_unpackhi_epi8(ymm_U_8dl, _mm256_setzero_si256());

					__m256i ymm_V_8dl = _mm256_unpacklo_epi8(ymm_V, ymm_V);
					__m256i ymm_V_8dh = _mm256_unpackhi_epi8(ymm_V, ymm_V);

					__m256i ymm_V_16l = _mm256_unpacklo_epi8(ymm_V_8dl, _mm256_setzero_si256());
					__m256i ymm_V_16h = _mm256_unpackhi_epi8(ymm_V_8dl, _mm256_setzero_si256());


					ymm_U0_32l = _mm256_unpacklo_epi16(ymm_U_16l, _mm256_setzero_si256());
					ymm_V0_32l = _mm256_unpacklo_epi16(ymm_V_16l, _mm256_setzero_si256());

					ymm_U0_32h = _mm256_unpackhi_epi16(ymm_U_16l, _mm256_setzero_si256());
					ymm_V0_32h = _mm256_unpackhi_epi16(ymm_V_16l, _mm256_setzero_si256());

					ymm_U1_32l = _mm256_unpacklo_epi16(ymm_U_16h, _mm256_setzero_si256());
					ymm_V1_32l = _mm256_unpacklo_epi16(ymm_V_16h, _mm256_setzero_si256());

					ymm_U1_32h = _mm256_unpackhi_epi16(ymm_U_16h, _mm256_setzero_si256());
					ymm_V1_32h = _mm256_unpackhi_epi16(ymm_V_16h, _mm256_setzero_si256());


				}
				/*				if (cs == 2)
								{
									if (!bCacheLoad)
									{
										_mm_prefetch((const CHAR*)(l_srcp_Y + 128), _MM_HINT_NTA);
										_mm_prefetch((const CHAR*)(l_srcp_Y + 128 + 64), _MM_HINT_NTA);

										_mm_prefetch((const CHAR*)(l_srcp_U + 64), _MM_HINT_NTA);
										_mm_prefetch((const CHAR*)(l_srcp_V + 64), _MM_HINT_NTA);
									}

									__m256i ymm_Y0_l = _mm256_load_si256((const __m256i*)l_srcp_Y); // should always load from 64-bit aligned start of row in AVS+ 3.7.3 (and later ?)
									__m256i ymm_Y0_h = _mm256_load_si256((const __m256i*)(l_srcp_Y + 32));
									__m256i ymm_Y1_l = _mm256_load_si256((const __m256i*)(l_srcp_Y + 64));
									__m256i ymm_Y1_h = _mm256_load_si256((const __m256i*)(l_srcp_Y + 96));

									__m256i ymm_U_l = _mm256_load_si256((const __m256i*)(l_srcp_U));
									__m256i ymm_U_h = _mm256_load_si256((const __m256i*)(l_srcp_U + 32));

									__m256i ymm_V_l = _mm256_load_si256((const __m256i*)(l_srcp_V));
									__m256i ymm_V_h = _mm256_load_si256((const __m256i*)(l_srcp_V + 32));

									ymm_Y0_16l = _mm256_permute2x128_si256(ymm_Y0_l, ymm_Y0_h, 0x20);
									ymm_Y0_16h = _mm256_permute2x128_si256(ymm_Y0_l, ymm_Y0_h, 0x31);

									ymm_Y1_16l = _mm256_permute2x128_si256(ymm_Y1_l, ymm_Y1_h, 0x20);
									ymm_Y1_16h = _mm256_permute2x128_si256(ymm_Y1_l, ymm_Y1_h, 0x31);

									if (bps == 10)
									{
										ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 2);
										ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 2);

										ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 2);
										ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 2);

										ymm_U_l = _mm256_srli_epi16(ymm_U_l, 2);
										ymm_U_h = _mm256_srli_epi16(ymm_U_h, 2);

										ymm_V_l = _mm256_srli_epi16(ymm_V_l, 2);
										ymm_V_h = _mm256_srli_epi16(ymm_V_h, 2);
									}

									if (bps == 12)
									{
										ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 4);
										ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 4);

										ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 4);
										ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 4);

										ymm_U_l = _mm256_srli_epi16(ymm_U_l, 4);
										ymm_U_h = _mm256_srli_epi16(ymm_U_h, 4);

										ymm_V_l = _mm256_srli_epi16(ymm_V_l, 4);
										ymm_V_h = _mm256_srli_epi16(ymm_V_h, 4);
									}

									if (bps == 14)
									{
										ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 6);
										ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 6);

										ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 6);
										ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 6);

										ymm_U_l = _mm256_srli_epi16(ymm_U_l, 6);
										ymm_U_h = _mm256_srli_epi16(ymm_U_h, 6);

										ymm_V_l = _mm256_srli_epi16(ymm_V_l, 6);
										ymm_V_h = _mm256_srli_epi16(ymm_V_h, 6);
									}

									if (bps == 16)
									{
										ymm_Y0_16l = _mm256_srli_epi16(ymm_Y0_16l, 8);
										ymm_Y1_16l = _mm256_srli_epi16(ymm_Y1_16l, 8);

										ymm_Y0_16h = _mm256_srli_epi16(ymm_Y0_16h, 8);
										ymm_Y1_16h = _mm256_srli_epi16(ymm_Y1_16h, 8);

										ymm_U_l = _mm256_srli_epi16(ymm_U_l, 8);
										ymm_U_h = _mm256_srli_epi16(ymm_U_h, 8);

										ymm_V_l = _mm256_srli_epi16(ymm_V_l, 8);
										ymm_V_h = _mm256_srli_epi16(ymm_V_h, 8);
									}


									ymm_U0_16l = _mm256_unpacklo_epi16(ymm_U_l, ymm_U_l);
									ymm_V0_16l = _mm256_unpacklo_epi16(ymm_V_l, ymm_V_l);

									ymm_U1_16l = _mm256_unpacklo_epi16(ymm_U_h, ymm_U_h);
									ymm_V1_16l = _mm256_unpacklo_epi16(ymm_V_h, ymm_V_h);

									ymm_U0_16h = _mm256_unpackhi_epi16(ymm_U_l, ymm_U_l);
									ymm_V0_16h = _mm256_unpackhi_epi16(ymm_V_l, ymm_V_l);

									ymm_U1_16h = _mm256_unpackhi_epi16(ymm_U_h, ymm_U_h);
									ymm_V1_16h = _mm256_unpackhi_epi16(ymm_V_h, ymm_V_h);


								}
								*/
				ymm_U0_32l = _mm256_sub_epi32(ymm_U0_32l, ymm_dw_cbias); // superscalarity of 3 at IceLake/SkyLake
				ymm_V0_32l = _mm256_sub_epi32(ymm_V0_32l, ymm_dw_cbias);

				ymm_U0_32h = _mm256_sub_epi32(ymm_U0_32h, ymm_dw_cbias);
				ymm_V0_32h = _mm256_sub_epi32(ymm_V0_32h, ymm_dw_cbias);

				ymm_U1_32l = _mm256_sub_epi32(ymm_U1_32l, ymm_dw_cbias); // superscalarity of 3 at IceLake/SkyLake
				ymm_V1_32l = _mm256_sub_epi32(ymm_V1_32l, ymm_dw_cbias);

				ymm_U1_32h = _mm256_sub_epi32(ymm_U1_32h, ymm_dw_cbias);
				ymm_V1_32h = _mm256_sub_epi32(ymm_V1_32h, ymm_dw_cbias);


				__m256i ymm_V_dl32l_addR = _mm256_mullo_epi32(ymm_V0_32l, ymm_dwKr);
				__m256i ymm_V_dl32h_addR = _mm256_mullo_epi32(ymm_V0_32h, ymm_dwKr);

				__m256i ymm_V_dh32l_addR = _mm256_mullo_epi32(ymm_V1_32l, ymm_dwKr);
				__m256i ymm_V_dh32h_addR = _mm256_mullo_epi32(ymm_V1_32h, ymm_dwKr);


				__m256i ymm_U_dl32l_addB = _mm256_mullo_epi32(ymm_U0_32l, ymm_dwKb);
				__m256i ymm_U_dl32h_addB = _mm256_mullo_epi32(ymm_U0_32h, ymm_dwKb);

				__m256i ymm_U_dh32l_addB = _mm256_mullo_epi32(ymm_U1_32l, ymm_dwKb);
				__m256i ymm_U_dh32h_addB = _mm256_mullo_epi32(ymm_U1_32h, ymm_dwKb);

				ymm_V_dl32l_addR = _mm256_srai_epi32(ymm_V_dl32l_addR, 13);
				ymm_V_dl32h_addR = _mm256_srai_epi32(ymm_V_dl32h_addR, 13);

				ymm_V_dh32l_addR = _mm256_srai_epi32(ymm_V_dh32l_addR, 13);
				ymm_V_dh32h_addR = _mm256_srai_epi32(ymm_V_dh32h_addR, 13);


				ymm_U_dl32l_addB = _mm256_srai_epi32(ymm_U_dl32l_addB, 13);
				ymm_U_dl32h_addB = _mm256_srai_epi32(ymm_U_dl32h_addB, 13);

				ymm_U_dh32l_addB = _mm256_srai_epi32(ymm_U_dh32l_addB, 13);
				ymm_U_dh32h_addB = _mm256_srai_epi32(ymm_U_dh32h_addB, 13);


				// R
				__m256i ymm_R_0_32l = _mm256_add_epi32(ymm_Y0_32l, ymm_V_dl32l_addR);
				__m256i ymm_R_0_32h = _mm256_add_epi32(ymm_Y0_32h, ymm_V_dl32h_addR);

				__m256i ymm_R_1_32l = _mm256_add_epi32(ymm_Y1_32l, ymm_V_dh32l_addR);
				__m256i ymm_R_1_32h = _mm256_add_epi32(ymm_Y1_32h, ymm_V_dh32h_addR);

				// B
				__m256i ymm_B_0_32l = _mm256_add_epi32(ymm_Y0_32l, ymm_U_dl32l_addB);
				__m256i ymm_B_0_32h = _mm256_add_epi32(ymm_Y0_32h, ymm_U_dl32h_addB);

				__m256i ymm_B_1_32l = _mm256_add_epi32(ymm_Y1_32l, ymm_U_dh32l_addB);
				__m256i ymm_B_1_32h = _mm256_add_epi32(ymm_Y1_32h, ymm_U_dh32h_addB);

				//G
				__m256i ymm_V_dl32l_subG = _mm256_mullo_epi32(ymm_V0_32l, ymm_dwKgv);
				__m256i ymm_V_dl32h_subG = _mm256_mullo_epi32(ymm_V0_32h, ymm_dwKgv);

				__m256i ymm_V_dh32l_subG = _mm256_mullo_epi32(ymm_V1_32l, ymm_dwKgv);
				__m256i ymm_V_dh32h_subG = _mm256_mullo_epi32(ymm_V1_32h, ymm_dwKgv);


				__m256i ymm_U_dl32l_subG = _mm256_mullo_epi32(ymm_U0_32l, ymm_dwKgu);
				__m256i ymm_U_dl32h_subG = _mm256_mullo_epi32(ymm_U0_32h, ymm_dwKgu);

				__m256i ymm_U_dh32l_subG = _mm256_mullo_epi32(ymm_U1_32l, ymm_dwKgu);
				__m256i ymm_U_dh32h_subG = _mm256_mullo_epi32(ymm_U1_32h, ymm_dwKgu);


				ymm_V_dl32l_subG = _mm256_srai_epi32(ymm_V_dl32l_subG, 13);
				ymm_V_dl32h_subG = _mm256_srai_epi32(ymm_V_dl32h_subG, 13);

				ymm_V_dh32l_subG = _mm256_srai_epi32(ymm_V_dh32l_subG, 13);
				ymm_V_dh32h_subG = _mm256_srai_epi32(ymm_V_dh32h_subG, 13);


				ymm_U_dl32l_subG = _mm256_srai_epi32(ymm_U_dl32l_subG, 13);
				ymm_U_dl32h_subG = _mm256_srai_epi32(ymm_U_dl32h_subG, 13);

				ymm_U_dh32l_subG = _mm256_srai_epi32(ymm_U_dh32l_subG, 13);
				ymm_U_dh32h_subG = _mm256_srai_epi32(ymm_U_dh32h_subG, 13);


				__m256i ymm_G_0_32l = _mm256_sub_epi32(ymm_Y0_32l, ymm_V_dl32l_subG);
				__m256i ymm_G_0_32h = _mm256_sub_epi32(ymm_Y0_32h, ymm_V_dl32h_subG);

				ymm_G_0_32l = _mm256_sub_epi32(ymm_G_0_32l, ymm_U_dl32l_subG);
				ymm_G_0_32h = _mm256_sub_epi32(ymm_G_0_32h, ymm_U_dl32h_subG);

				__m256i ymm_G_1_32l = _mm256_sub_epi32(ymm_Y1_32l, ymm_V_dh32l_subG);
				__m256i ymm_G_1_32h = _mm256_sub_epi32(ymm_Y1_32h, ymm_V_dh32h_subG);

				ymm_G_1_32l = _mm256_sub_epi32(ymm_G_1_32l, ymm_U_dh32l_subG);
				ymm_G_1_32h = _mm256_sub_epi32(ymm_G_1_32h, ymm_U_dh32h_subG);

				// RGB post processing with gain and offset 
				__m256i ymm_RGBoffset = _mm256_set1_epi32(RGBo);
				__m256i ymm_RGBgain = _mm256_set1_epi32(RGBg);

				ymm_R_0_32l = _mm256_add_epi32(ymm_R_0_32l, ymm_RGBoffset);
				ymm_R_0_32h = _mm256_add_epi32(ymm_R_0_32h, ymm_RGBoffset);
				ymm_R_1_32l = _mm256_add_epi32(ymm_R_1_32l, ymm_RGBoffset);
				ymm_R_1_32h = _mm256_add_epi32(ymm_R_1_32h, ymm_RGBoffset);

				ymm_G_0_32l = _mm256_add_epi32(ymm_G_0_32l, ymm_RGBoffset);
				ymm_G_0_32h = _mm256_add_epi32(ymm_G_0_32h, ymm_RGBoffset);
				ymm_G_1_32l = _mm256_add_epi32(ymm_G_1_32l, ymm_RGBoffset);
				ymm_G_1_32h = _mm256_add_epi32(ymm_G_1_32h, ymm_RGBoffset);

				ymm_B_0_32l = _mm256_add_epi32(ymm_B_0_32l, ymm_RGBoffset);
				ymm_B_0_32h = _mm256_add_epi32(ymm_B_0_32h, ymm_RGBoffset);
				ymm_B_1_32l = _mm256_add_epi32(ymm_B_1_32l, ymm_RGBoffset);
				ymm_B_1_32h = _mm256_add_epi32(ymm_B_1_32h, ymm_RGBoffset);


				ymm_R_0_32l = _mm256_mullo_epi32(ymm_R_0_32l, ymm_RGBgain);
				ymm_R_0_32h = _mm256_mullo_epi32(ymm_R_0_32h, ymm_RGBgain);
				ymm_R_1_32l = _mm256_mullo_epi32(ymm_R_1_32l, ymm_RGBgain);
				ymm_R_1_32h = _mm256_mullo_epi32(ymm_R_1_32h, ymm_RGBgain);

				ymm_G_0_32l = _mm256_mullo_epi32(ymm_G_0_32l, ymm_RGBgain);
				ymm_G_0_32h = _mm256_mullo_epi32(ymm_G_0_32h, ymm_RGBgain);
				ymm_G_1_32l = _mm256_mullo_epi32(ymm_G_1_32l, ymm_RGBgain);
				ymm_G_1_32h = _mm256_mullo_epi32(ymm_G_1_32h, ymm_RGBgain);

				ymm_B_0_32l = _mm256_mullo_epi32(ymm_B_0_32l, ymm_RGBgain);
				ymm_B_0_32h = _mm256_mullo_epi32(ymm_B_0_32h, ymm_RGBgain);
				ymm_B_1_32l = _mm256_mullo_epi32(ymm_B_1_32l, ymm_RGBgain);
				ymm_B_1_32h = _mm256_mullo_epi32(ymm_B_1_32h, ymm_RGBgain);


				ymm_R_0_32l = _mm256_srai_epi32(ymm_R_0_32l, RGB_DIV_SHIFT_32);
				ymm_R_0_32h = _mm256_srai_epi32(ymm_R_0_32h, RGB_DIV_SHIFT_32);
				ymm_R_1_32l = _mm256_srai_epi32(ymm_R_1_32l, RGB_DIV_SHIFT_32);
				ymm_R_1_32h = _mm256_srai_epi32(ymm_R_1_32h, RGB_DIV_SHIFT_32);

				ymm_G_0_32l = _mm256_srai_epi32(ymm_G_0_32l, RGB_DIV_SHIFT_32);
				ymm_G_0_32h = _mm256_srai_epi32(ymm_G_0_32h, RGB_DIV_SHIFT_32);
				ymm_G_1_32l = _mm256_srai_epi32(ymm_G_1_32l, RGB_DIV_SHIFT_32);
				ymm_G_1_32h = _mm256_srai_epi32(ymm_G_1_32h, RGB_DIV_SHIFT_32);

				ymm_B_0_32l = _mm256_srai_epi32(ymm_B_0_32l, RGB_DIV_SHIFT_32);
				ymm_B_0_32h = _mm256_srai_epi32(ymm_B_0_32h, RGB_DIV_SHIFT_32);
				ymm_B_1_32l = _mm256_srai_epi32(ymm_B_1_32l, RGB_DIV_SHIFT_32);
				ymm_B_1_32h = _mm256_srai_epi32(ymm_B_1_32h, RGB_DIV_SHIFT_32);


				//pack 32bit to 8bit
				__m256i ymm_R_0_16 = _mm256_packus_epi32(ymm_R_0_32l, ymm_R_0_32h);
				__m256i ymm_R_1_16 = _mm256_packus_epi32(ymm_R_1_32l, ymm_R_1_32h);

				__m256i ymm_R_8 = _mm256_packus_epi16(ymm_R_0_16, ymm_R_1_16);

				__m256i ymm_G_0_16 = _mm256_packus_epi32(ymm_G_0_32l, ymm_G_0_32h);
				__m256i ymm_G_1_16 = _mm256_packus_epi32(ymm_G_1_32l, ymm_G_1_32h);

				__m256i ymm_G_8 = _mm256_packus_epi16(ymm_G_0_16, ymm_G_1_16);

				__m256i ymm_B_0_16 = _mm256_packus_epi32(ymm_B_0_32l, ymm_B_0_32h);
				__m256i ymm_B_1_16 = _mm256_packus_epi32(ymm_B_1_32l, ymm_B_1_32h);

				__m256i ymm_B_8 = _mm256_packus_epi16(ymm_B_0_16, ymm_B_1_16);

				/*
				if (bCacheStore)
				{
					_mm256_store_si256((__m256i*)l_dstp_R, ymm_R_8);
					_mm256_store_si256((__m256i*)l_dstp_G, ymm_G_8);
					_mm256_store_si256((__m256i*)l_dstp_B, ymm_B_8);
				}
				else
				{
					_mm256_stream_si256((__m256i*)l_dstp_R, ymm_R_8);
					_mm256_stream_si256((__m256i*)l_dstp_G, ymm_G_8);
					_mm256_stream_si256((__m256i*)l_dstp_B, ymm_B_8);

				}
				*/

				// convert planar RGB to BGRA32

				__m256i ymm_blend_G = _mm256_setr_epi8(0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0);
				__m256i ymm_blend_R = _mm256_setr_epi8(0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0);

				__m256i ymm_shuf_B_0_7 = _mm256_setr_epi8(0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255, 255);
				__m256i ymm_shuf_B_8_15 = _mm256_setr_epi8(8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255, 255);

				__m256i ymm_shuf_G_0_7 = _mm256_setr_epi8(255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255);
				__m256i ymm_shuf_G_8_15 = _mm256_setr_epi8(255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255);

				__m256i ymm_shuf_R_0_7 = _mm256_setr_epi8(255, 255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255);
				__m256i ymm_shuf_R_8_15 = _mm256_setr_epi8(255, 255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255);


				__m256i ymm_R_0_15 = _mm256_permute4x64_epi64(ymm_R_8, 0x44);
				__m256i ymm_G_0_15 = _mm256_permute4x64_epi64(ymm_G_8, 0x44);
				__m256i ymm_B_0_15 = _mm256_permute4x64_epi64(ymm_B_8, 0x44);

				// B
				__m256i ymm_BGRA_0_7 = _mm256_shuffle_epi8(ymm_B_0_15, ymm_shuf_B_0_7);
				__m256i ymm_BGRA_8_15 = _mm256_shuffle_epi8(ymm_B_0_15, ymm_shuf_B_8_15);

				// G
				__m256i ymm_G_0_7 = _mm256_shuffle_epi8(ymm_G_0_15, ymm_shuf_G_0_7);
				__m256i ymm_G_8_15 = _mm256_shuffle_epi8(ymm_G_0_15, ymm_shuf_G_8_15);

				ymm_BGRA_0_7 = _mm256_blendv_epi8(ymm_BGRA_0_7, ymm_G_0_7, ymm_blend_G);
				ymm_BGRA_8_15 = _mm256_blendv_epi8(ymm_BGRA_8_15, ymm_G_8_15, ymm_blend_G);

				//R
				__m256i ymm_R_0_7 = _mm256_shuffle_epi8(ymm_R_0_15, ymm_shuf_R_0_7);
				__m256i ymm_R_8_15 = _mm256_shuffle_epi8(ymm_R_0_15, ymm_shuf_R_8_15);

				ymm_BGRA_0_7 = _mm256_blendv_epi8(ymm_BGRA_0_7, ymm_R_0_7, ymm_blend_R);
				ymm_BGRA_8_15 = _mm256_blendv_epi8(ymm_BGRA_8_15, ymm_R_8_15, ymm_blend_R);

				////// BGR 16_23 and 24_31
				__m256i ymm_R_16_31 = _mm256_permute4x64_epi64(ymm_R_8, 0xEE);
				__m256i ymm_G_16_31 = _mm256_permute4x64_epi64(ymm_G_8, 0xEE);
				__m256i ymm_B_16_31 = _mm256_permute4x64_epi64(ymm_B_8, 0xEE);

				// B
				__m256i ymm_BGRA_16_23 = _mm256_shuffle_epi8(ymm_B_16_31, ymm_shuf_B_0_7);
				__m256i ymm_BGRA_24_31 = _mm256_shuffle_epi8(ymm_B_16_31, ymm_shuf_B_8_15);

				// G
				__m256i ymm_G_16_23 = _mm256_shuffle_epi8(ymm_G_16_31, ymm_shuf_G_0_7);
				__m256i ymm_G_24_31 = _mm256_shuffle_epi8(ymm_G_16_31, ymm_shuf_G_8_15);

				ymm_BGRA_16_23 = _mm256_blendv_epi8(ymm_BGRA_16_23, ymm_G_16_23, ymm_blend_G);
				ymm_BGRA_24_31 = _mm256_blendv_epi8(ymm_BGRA_24_31, ymm_G_24_31, ymm_blend_G);

				//R
				__m256i ymm_R_16_23 = _mm256_shuffle_epi8(ymm_R_16_31, ymm_shuf_R_0_7);
				__m256i ymm_R_24_31 = _mm256_shuffle_epi8(ymm_R_16_31, ymm_shuf_R_8_15);

				ymm_BGRA_16_23 = _mm256_blendv_epi8(ymm_BGRA_16_23, ymm_R_16_23, ymm_blend_R);
				ymm_BGRA_24_31 = _mm256_blendv_epi8(ymm_BGRA_24_31, ymm_R_24_31, ymm_blend_R);

				if (bCacheStore)
				{
					_mm256_store_si256((__m256i*)l_dstp_BGRA, ymm_BGRA_0_7);
					_mm256_store_si256((__m256i*)(l_dstp_BGRA + 32), ymm_BGRA_8_15);
					_mm256_store_si256((__m256i*)(l_dstp_BGRA + 64), ymm_BGRA_16_23);
					_mm256_store_si256((__m256i*)(l_dstp_BGRA + 96), ymm_BGRA_24_31);
				}
				else
				{
					_mm256_stream_si256((__m256i*)l_dstp_BGRA, ymm_BGRA_0_7);
					_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 32), ymm_BGRA_8_15);
					_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 64), ymm_BGRA_16_23);
					_mm256_stream_si256((__m256i*)(l_dstp_BGRA + 96), ymm_BGRA_24_31);
				}

				if (cs == 1)
				{
					l_srcp_Y += 32; // in bytes
					l_srcp_U += 16;
					l_srcp_V += 16;
				}
				else
				{
					l_srcp_Y += 64; // in bytes
					l_srcp_U += 32;
					l_srcp_V += 32;
				}

				l_dstp_R += 32;
				l_dstp_G += 32;
				l_dstp_B += 32;

				l_dstp_BGRA += 32 * 4;
			}

			// last cols
			int iUVadv = 0;
			for (int col = col32; col < row_size_Ysamples; ++col)
			{
				int iY;
				int iU;
				int iV;

				if (cs == 1)
				{
					iY = *l_srcp_Y;
					iU = *l_srcp_U - 128;
					iV = *l_srcp_V - 128;
				}
				else
				{
					iY = (int) * ((unsigned short*)l_srcp_Y);
					iU = (int) * ((unsigned short*)l_srcp_U);
					iV = (int) * ((unsigned short*)l_srcp_V);

					if (bps == 10)
					{
						iY = iY >> 2;
						iU = iU >> 2;
						iV = iV >> 2;
					}

					if (bps == 12)
					{
						iY = iY >> 4;
						iU = iU >> 4;
						iV = iV >> 4;
					}

					if (bps == 14)
					{
						iY = iY >> 6;
						iU = iU >> 6;
						iV = iV >> 6;
					}

					if (bps == 16)
					{
						iY = iY >> 8;
						iU = iU >> 8;
						iV = iV >> 8;
					}

					iU = iU - 128;
					iV = iV - 128;

				}

				int iR = iY + ((iV * Kr) >> 13);
				int iB = iY + ((iU * Kb) >> 13);
				int iG = iY - ((iU * Kgu) >> 13) - ((iV * Kgv) >> 13);

				iR = ((iR + RGBo) * RGBg) >> RGB_DIV_SHIFT_32;
				iG = ((iG + RGBo) * RGBg) >> RGB_DIV_SHIFT_32;
				iB = ((iB + RGBo) * RGBg) >> RGB_DIV_SHIFT_32;

				if (iR < 0) iR = 0; if (iR > 255) iR = 255;
				if (iG < 0) iG = 0; if (iG > 255) iG = 255;
				if (iB < 0) iB = 0; if (iB > 255) iB = 255;

				int iBGRA = 0 | (unsigned char)iR << 16 | (unsigned char)iG << 8 | (unsigned char)iB;
				*(int*)(l_dstp_BGRA) = iBGRA;


				if (cs == 1)
				{
					l_srcp_Y += 1; // in bytes
					if (iUVadv % 2 != 0)
					{
						l_srcp_U += 1;
						l_srcp_V += 1;
					}
				}
				else
				{
					l_srcp_Y += 2; // in bytes
					if (iUVadv % 2 != 0)
					{
						l_srcp_U += 2;
						l_srcp_V += 2;
					}
				}

				iUVadv++;
				l_dstp_BGRA += 4;

			}

		}


	}

	if (!bCacheStore)
		_mm_sfence();
}


AVSValue __cdecl Create_Decode(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	return new DecodeYUVtoRGB(args[0].AsClip(), args[1].AsInt(1), args[2].AsInt(0), args[3].AsInt(64), args[4].AsInt(0), args[5].AsBool(true), args[6].AsBool(false), args[7].AsInt(16), env);
}

const AVS_Linkage* AVS_linkage = 0;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
	AVS_linkage = vectors;
	env->AddFunction("DecodeYUVtoRGB", "c[threads]i[matrix]i[gain]i[offset]i[cl]b[cs]b[ib]i", Create_Decode, 0);

	return "Decode YUV to RGB sample plugin";
}


