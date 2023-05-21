#include "include\avisynth.h"
#include <windows.h>
#include <immintrin.h>

#define RGB_DIV_SHIFT 6

class DecodeYUVtoRGB : public GenericVideoFilter
{
	int threads;
	int _cpuFlags;

	short Kr; // div64
	short Kb; // div64
	short Kgu; // div64
	short Kgv; // div64

	int Matrix;

	short RGBgain;
	short RGBoffset;

	bool bCacheLoad;
	bool bCacheStore;

public:
	DecodeYUVtoRGB(PClip _child, int threads_, int _matrix, int _gain, int _offset, bool _cl, bool _cs, IScriptEnvironment* env) : GenericVideoFilter(_child),
		threads(threads_),
		Matrix(_matrix),
		RGBgain((short)_gain),
		RGBoffset((short)_offset),
		bCacheLoad(_cl),
		bCacheStore(_cs)
	{
		_cpuFlags = env->GetCPUFlags();
		//		vi.pixel_type = VideoInfo::CS_RGBP8;
		vi.pixel_type = VideoInfo::CS_BGR32;

		if (!(_cpuFlags & CPUF_AVX2))
		{
			env->ThrowError("DecodeYUVtoRGB: Only AVX2 and later SIMD co-processor supported.");
		}

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

	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env)
	{
		PVideoFrame dst = env->NewVideoFrame(vi);
		VideoInfo vi_src = child->GetVideoInfo();
		PVideoFrame src = child->GetFrame(n, env);

		DWORD dwOldProtY;
		DWORD dwOldProtU;
		DWORD dwOldProtV;
		bool bRes;
		
		if (vi_src.ComponentSize() == 1)
		{

			if (bCacheLoad && bCacheStore)
				DecodeYUV420<true, true, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && bCacheStore)
				DecodeYUV420<false, true, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (bCacheLoad && !bCacheStore)
				DecodeYUV420<true, false, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && !bCacheStore)
				DecodeYUV420<false, false, 8, 1>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

		}
		else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 10))
		{
			if (bCacheLoad && bCacheStore)
				DecodeYUV420<true, true, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && bCacheStore)
				DecodeYUV420<false, true, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (bCacheLoad && !bCacheStore)
				DecodeYUV420<true, false, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && !bCacheStore)
				DecodeYUV420<false, false, 10, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
		}
		else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 12))
		{
			if (bCacheLoad && bCacheStore)
				DecodeYUV420<true, true, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && bCacheStore)
				DecodeYUV420<false, true, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (bCacheLoad && !bCacheStore)
				DecodeYUV420<true, false, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && !bCacheStore)
				DecodeYUV420<false, false, 12, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
		}
		else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 14))
		{
			if (bCacheLoad && bCacheStore)
				DecodeYUV420<true, true, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && bCacheStore)
				DecodeYUV420<false, true, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (bCacheLoad && !bCacheStore)
				DecodeYUV420<true, false, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && !bCacheStore)
				DecodeYUV420<false, false, 14, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
		}
		else if ((vi_src.ComponentSize() == 2) && (vi_src.BitsPerComponent() == 16))
		{
			if (bCacheLoad && bCacheStore)
				DecodeYUV420<true, true, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && bCacheStore)
				DecodeYUV420<false, true, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (bCacheLoad && !bCacheStore)
				DecodeYUV420<true, false, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);

			if (!bCacheLoad && !bCacheStore)
				DecodeYUV420<false, false, 16, 2>(dst, src, vi, vi_src, Kr, Kb, Kgu, Kgv, RGBgain, RGBoffset, threads, _cpuFlags);
		}
		else
			env->ThrowError("DecodeYV12toRGB: Only 8bit and 10bit input supported.");

		return dst;
	}

	template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
	void DecodeYUV420(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags);

};


template <bool bCacheLoad, bool bCacheStore, int bps, int cs>
void DecodeYUVtoRGB::DecodeYUV420(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, short Kr, short Kb, short Kgu, short Kgv, short RGBg, short RGBo, int threads, int cpuFlags)
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
				const __m256i ymm_w_cbias = _mm256_set1_epi16(128);

				const __m256i ymm_wKr = _mm256_set1_epi16(Kr); 
				const __m256i ymm_wKb = _mm256_set1_epi16(Kb); 

				const __m256i ymm_wKgu = _mm256_set1_epi16(Kgu); 
				const __m256i ymm_wKgv = _mm256_set1_epi16(Kgv);

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
					const __m256i ymm_RGBoffset = _mm256_set1_epi16(RGBo);
					const __m256i ymm_RGBgain = _mm256_set1_epi16(RGBg);

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

					const __m256i ymm_blend_G = _mm256_setr_epi8(0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0);
					const __m256i ymm_blend_R = _mm256_setr_epi8(0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0);

					const __m256i ymm_shuf_B_0_7 = _mm256_setr_epi8(0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255, 255);
					const __m256i ymm_shuf_B_8_15 = _mm256_setr_epi8(8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255, 255);

					const __m256i ymm_shuf_G_0_7 = _mm256_setr_epi8(255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255, 255);
					const __m256i ymm_shuf_G_8_15 = _mm256_setr_epi8(255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255, 255);

					const __m256i ymm_shuf_R_0_7 = _mm256_setr_epi8(255, 255, 0, 255, 255, 255, 1, 255, 255, 255, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255, 5, 255, 255, 255, 6, 255, 255, 255, 7, 255);
					const __m256i ymm_shuf_R_8_15 = _mm256_setr_epi8(255, 255, 8, 255, 255, 255, 9, 255, 255, 255, 10, 255, 255, 255, 11, 255, 255, 255, 12, 255, 255, 255, 13, 255, 255, 255, 14, 255, 255, 255, 15, 255);


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

					int iR = iY + ((iV * Kr) >> 6);
					int iB = iY + ((iU * Kb) >> 6);
					int iG = iY - ((iU * Kgu) >> 6) - ((iV * Kgv) >> 6);

					iR = ((iR + RGBo) * RGBg) >> RGB_DIV_SHIFT;
					iG = ((iG + RGBo) * RGBg) >> RGB_DIV_SHIFT;
					iB = ((iB + RGBo) * RGBg) >> RGB_DIV_SHIFT;

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
	return new DecodeYUVtoRGB(args[0].AsClip(), args[1].AsInt(1), args[2].AsInt(0), args[3].AsInt(64), args[4].AsInt(0), args[5].AsBool(true), args[6].AsBool(false), env);
}

const AVS_Linkage* AVS_linkage = 0;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
	AVS_linkage = vectors;
	env->AddFunction("DecodeYUVtoRGB", "c[threads]i[matrix]i[gain]i[offset]i[cl]b[cs]b", Create_Decode, 0);

	return "Decode YUV to RGB sample plugin";
}


