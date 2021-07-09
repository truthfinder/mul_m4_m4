//#include <windows.h>
//#include <conio.h>

#ifdef IACA_MARKS_OFF

#	define IACA_FROM
#	define IACA_TO

#else

#if defined IACA3
#	include "../iaca/iacaMarks.h"
#elif defined IACA2
#	include "../iaca2/iacaMarks.h"
#else
#	error iaca version not defined
#endif

#	ifdef _WIN64
#		define IACA_FROM IACA_VC64_START
#		define IACA_TO IACA_VC64_END
#	else
#		define IACA_FROM IACA_START
#		define IACA_TO IACA_END
#	endif

#endif

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <chrono>
#include <array>
#include <stdlib.h>

#define SKIP_GLM 1

#if !defined SKIP_GLM
#define GLM_FORCE_SSE2
#define GLM_FORCE_SWIZZLE
#include <glm/vec3.hpp>   // glm::vec3
#include <glm/vec4.hpp>   // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/simd/matrix.h>
#endif

#ifdef _MSC_VER
#include <intrin.h>
#define __vmathcall  __vectorcall
#else
#include <cpuid.h>
#include <x86intrin.h>
#define __vmathcall
#endif

//std::chrono::duration_cast<std::chrono::seconds>(tp-tp).count()

template <typename T> class check_ts { // temporary class for performance debugging purposes
public:
	static int64_t freq() noexcept { return T::freq(); }

	check_ts(bool const started = true) noexcept : t(0ll), r(0ll), f(1. / T::freq()) { if (started) { start(); } }
	~check_ts() noexcept { stop(); }

	int64_t      start() noexcept { return t = T::tics_start(); }
	int64_t       stop() noexcept { return r = T::tics_stop() - t; }
	double    stop_sec() noexcept { return f * stop(); }
	int64_t      reset() noexcept { int64_t c = T::tics_start(); r = c - t; t = c; return r; }
	double   reset_sec() noexcept { return f * reset(); }
	int64_t     passed() const noexcept { return r; }
	int64_t  passed_ms() const noexcept {
		std::chrono::duration<int64_t, std::ratio<1, T::freq()>> const u(passed());
		return std::chrono::duration_cast<std::chrono::milliseconds>(u).count();
	}
	int64_t passed_mks() const noexcept {
		std::chrono::duration<int64_t, std::ratio<1, T::freq()>> const u(passed());
		return std::chrono::duration_cast<std::chrono::microseconds>(u).count();
	}
	int64_t passed_ns() const noexcept {
		return passed();
	}
	double  passed_sec() const noexcept { return f * passed(); }
	int64_t    current() const noexcept { return T::tics() - t; }
	double current_sec() const noexcept { return f * current(); }

	check_ts& operator += (check_ts const& ts) noexcept { r += ts.r; return *this; }
	check_ts& operator -= (check_ts const& ts) noexcept { r -= ts.r; return *this; }
	check_ts& operator += (int64_t const   d) noexcept { r += d; return *this; }
	check_ts& operator -= (int64_t const   d) noexcept { r -= d; return *this; }

private:
	int64_t t, r;
	double f;
};

class hrc_timer { // temporary class for performance debugging purposes
	using clock = std::chrono::high_resolution_clock;
	//using frq = std::chrono::nanoseconds = clock::type?;
public:
	static int64_t tics_start() noexcept { return clock::now().time_since_epoch().count(); }
	static int64_t tics_stop() noexcept { return tics_start(); }
	static constexpr int64_t freq() noexcept { return clock::period::den; }
};

class rdtsc_timer {
public:
	static int64_t tics_start() noexcept {
		int info[4];
#		ifdef WIN32
		__cpuid(info, 1);
#		else
		__cpuid(1, info[0], info[1], info[2], info[3]);
#		endif
		return __rdtsc();
	}
	static int64_t tics_stop() noexcept {
		int info[4];
		unsigned int aux;
		int64_t r = __rdtscp(&aux);
#		ifdef WIN32
		__cpuid(info, 1);
#		else
		__cpuid(1, info[0], info[1], info[2], info[3]);
#		endif
		return r;
	}
	static int64_t freq() noexcept {
		return 3'700'000'000ll;
	}
};

typedef check_ts<hrc_timer> hts;
typedef check_ts<rdtsc_timer> rts;

//#define __FMA__
//#define __AVX2__
//#define __AVX512F__

// 00, 01, 02, 03
// 10, 11, 12, 13
// 20, 21, 22, 23
// 30, 31, 32, 33

// 00, 10, 20, 30
// 01, 11, 21, 31
// 02, 12, 22, 32
// 03, 13, 23, 33

struct alignas(sizeof(__m128)) vec4 {
	union {
		struct { float w, z, y, x; };
		__m128 fmm;
		float arr[4];
	};

	vec4() = default;
	vec4(float a, float b, float c, float d) : w(d), z(c), y(b), x(a) {}

	static bool equ(float const a, float const b, float t = .00001f) {
		return fabs(a-b) < t;
	}

	bool operator == (vec4 const& v) const {
		return equ(x, v.x) && equ(y, v.y) && equ(z, v.z) && equ(w, v.w);
	}
};

std::ostream& operator << (std::ostream& os, vec4 const& v) {
	return os
		<< std::setiosflags(std::ios_base::showpoint|std::ios_base::fixed)
		<< std::setprecision(3)
		<< '(' << v.x << ',' << v.y << ',' << v.z << ',' << v.w << ')';
}

struct alignas(sizeof(__m512)) mtx4 {
	union {
		struct {
			float
				_30, _20, _10, _00,
				_31, _21, _11, _01,
				_32, _22, _12, _02,
				_33, _23, _13, _03;
		};
		__m128 xmm[4];
		__m256 ymm[2];
		__m512 zmm;
		vec4 v[4];
	};

	mtx4() = default;
	mtx4(
		float i00, float i01, float i02, float i03,
		float i10, float i11, float i12, float i13,
		float i20, float i21, float i22, float i23,
		float i30, float i31, float i32, float i33)
		: _00(i00),  _01(i01),  _02(i02),  _03(i03)
		, _10(i10),  _11(i11),  _12(i12),  _13(i13)
		, _20(i20),  _21(i21),  _22(i22),  _23(i23)
		, _30(i30),  _31(i31),  _32(i32),  _33(i33)
	{}

	operator __m128 const* () const { return xmm; }
	operator __m128* () { return xmm; }

	bool operator == (mtx4 const& m) const {
		return v[0]==m.v[0] && v[1]==m.v[1] && v[2]==m.v[2] && v[3]==m.v[3];
	}

	static mtx4 identity() {
		return mtx4(
			1.0f, 0.f, 0.f, 0.f,
			0.0f, 1.f, 0.f, 0.f,
			0.0f, 0.f, 1.f, 0.f,
			0.0f, 0.f, 0.f, 1.f);
	}

	static mtx4 zero() {
		return mtx4(
			0.0f, 0.f, 0.f, 0.f,
			0.0f, 0.f, 0.f, 0.f,
			0.0f, 0.f, 0.f, 0.f,
			0.0f, 0.f, 0.f, 0.f);
	}
};

std::ostream& operator << (std::ostream& os, mtx4 const& m) {
	static const int w = 7;
	return os
		<< std::setiosflags(std::ios_base::showpoint|std::ios_base::fixed)
		<< std::setprecision(3)
		<< '|' << std::setw(w) << m._00 << ',' << std::setw(w) << m._01 << ',' << std::setw(w) << m._02 << ',' << std::setw(w) << m._03 << '|' << std::endl
		<< '|' << std::setw(w) << m._10 << ',' << std::setw(w) << m._11 << ',' << std::setw(w) << m._12 << ',' << std::setw(w) << m._13 << '|' << std::endl
		<< '|' << std::setw(w) << m._20 << ',' << std::setw(w) << m._21 << ',' << std::setw(w) << m._22 << ',' << std::setw(w) << m._23 << '|' << std::endl
		<< '|' << std::setw(w) << m._30 << ',' << std::setw(w) << m._31 << ',' << std::setw(w) << m._32 << ',' << std::setw(w) << m._33 << '|';
}

#if defined GLM_VERSION

std::ostream& operator << (std::ostream& os, glm::mat4 const& m) {
	return os
		<< '|' << m[0][0] << ',' << m[0][1] << ',' << m[0][2] << ',' << m[0][3] << '|' << std::endl
		<< '|' << m[1][0] << ',' << m[1][1] << ',' << m[1][2] << ',' << m[1][3] << '|' << std::endl
		<< '|' << m[2][0] << ',' << m[2][1] << ',' << m[2][2] << ',' << m[2][3] << '|' << std::endl
		<< '|' << m[3][0] << ',' << m[3][1] << ',' << m[3][2] << ',' << m[3][3] << '|';
}

void __vmathcall mul_mtx4_mtx4_glm(__m128* const _r, __m128 const* const _m, __m128 const* const _n) {
	static_assert(sizeof(glm::mat4) == sizeof(__m128) * 4);
	glm::mat4& r = **reinterpret_cast<glm::mat4* const*>(&_r);
	glm::mat4 const& m = **reinterpret_cast<glm::mat4 const* const*>(&_m);
	glm::mat4 const& n = **reinterpret_cast<glm::mat4 const* const*>(&_n);
	r = m * n;
}

void __vmathcall mul_mtx4_mtx4_glm_sse(__m128* const _r, __m128 const* const _m, __m128 const* const _n) {
	glm_mat4_mul(_m, _n, _r);
}

#endif // SKIP_GLM

void __vmathcall mul_mtx4_mtx4_loop(__m128* const _r, __m128 const* const _m, __m128 const* const _n) {
	typedef float raw4x4[4][4];

	raw4x4 const& m = **reinterpret_cast<raw4x4 const* const*>(&_m);
	raw4x4 const& n = **reinterpret_cast<raw4x4 const* const*>(&_n);
	raw4x4&       r = **reinterpret_cast<raw4x4* const*>(&_r);

	// _30, _20, _10, _00,
	// _31, _21, _11, _01,
	// _32, _22, _12, _02,
	// _33, _23, _13, _03;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			//r[i][j] = 0.f; // reference algorithm
			r[j][3-i] = 0.f;
			for (int k = 0; k < 4; ++k) {
				//IACA_FROM
				// 3.65*64~256; 2.97*64~192
				//r[i][j] += m[i][k] * n[k][j]; // reference algorithm
				//r[j][i] += m[k][i] * n[j][k]; // swap row/col due to transpose
				r[j][3-i] += m[k][3-i] * n[j][3-k]; // 3-col due to swapped row data
			}
			//IACA_TO
		}
	}
}

void __vmathcall mul_mtx4_mtx4_unroll(__m128* const _r, __m128 const* const _m, __m128 const* const _n) {
	mtx4 const& m = **reinterpret_cast<mtx4 const* const*>(&_m);
	mtx4 const& n = **reinterpret_cast<mtx4 const* const*>(&_n);
	mtx4&       r = **reinterpret_cast<mtx4* const*>(&_r);

	//IACA_FROM
	// 69.95, 64

	r._00 = m._00*n._00 + m._01*n._10 + m._02*n._20 + m._03*n._30;
	r._01 = m._00*n._01 + m._01*n._11 + m._02*n._21 + m._03*n._31;
	r._02 = m._00*n._02 + m._01*n._12 + m._02*n._22 + m._03*n._32;
	r._03 = m._00*n._03 + m._01*n._13 + m._02*n._23 + m._03*n._33;

	r._10 = m._10*n._00 + m._11*n._10 + m._12*n._20 + m._13*n._30;
	r._11 = m._10*n._01 + m._11*n._11 + m._12*n._21 + m._13*n._31;
	r._12 = m._10*n._02 + m._11*n._12 + m._12*n._22 + m._13*n._32;
	r._13 = m._10*n._03 + m._11*n._13 + m._12*n._23 + m._13*n._33;

	r._20 = m._20*n._00 + m._21*n._10 + m._22*n._20 + m._23*n._30;
	r._21 = m._20*n._01 + m._21*n._11 + m._22*n._21 + m._23*n._31;
	r._22 = m._20*n._02 + m._21*n._12 + m._22*n._22 + m._23*n._32;
	r._23 = m._20*n._03 + m._21*n._13 + m._22*n._23 + m._23*n._33;

	r._30 = m._30*n._00 + m._31*n._10 + m._32*n._20 + m._33*n._30;
	r._31 = m._30*n._01 + m._31*n._11 + m._32*n._21 + m._33*n._31;
	r._32 = m._30*n._02 + m._31*n._12 + m._32*n._22 + m._33*n._32;
	r._33 = m._30*n._03 + m._31*n._13 + m._32*n._23 + m._33*n._33;

	//IACA_TO
}

void __vmathcall mul_mtx4_mtx4_sse_v1(__m128* const r, __m128 const* const m, __m128 const* const n) {
	// 00, 01, 02, 03
	// 10, 11, 12, 13
	// 20, 21, 22, 23
	// 30, 31, 32, 33

	// 00, 10, 20, 30
	// 01, 11, 21, 31
	// 02, 12, 22, 32
	// 03, 13, 23, 33

	/* переставляем строчки текста
	r._00 = m._00*n._00 + m._01*n._10 + m._02*n._20 + m._03*n._30;
	r._10 = m._10*n._00 + m._11*n._10 + m._12*n._20 + m._13*n._30;
	r._20 = m._20*n._00 + m._21*n._10 + m._22*n._20 + m._23*n._30;
	r._30 = m._30*n._00 + m._31*n._10 + m._32*n._20 + m._33*n._30;

	r._01 = m._00*n._01 + m._01*n._11 + m._02*n._21 + m._03*n._31;
	r._11 = m._10*n._01 + m._11*n._11 + m._12*n._21 + m._13*n._31;
	r._21 = m._20*n._01 + m._21*n._11 + m._22*n._21 + m._23*n._31;
	r._31 = m._30*n._01 + m._31*n._11 + m._32*n._21 + m._33*n._31;

	r._02 = m._00*n._02 + m._01*n._12 + m._02*n._22 + m._03*n._32;
	r._12 = m._10*n._02 + m._11*n._12 + m._12*n._22 + m._13*n._32;
	r._22 = m._20*n._02 + m._21*n._12 + m._22*n._22 + m._23*n._32;
	r._32 = m._30*n._02 + m._31*n._12 + m._32*n._22 + m._33*n._32;

	r._03 = m._00*n._03 + m._01*n._13 + m._02*n._23 + m._03*n._33;
	r._13 = m._10*n._03 + m._11*n._13 + m._12*n._23 + m._13*n._33;
	r._23 = m._20*n._03 + m._21*n._13 + m._22*n._23 + m._23*n._33;
	r._33 = m._30*n._03 + m._31*n._13 + m._32*n._23 + m._33*n._33;
	*/

	// m[0] = {m._00, m._10, m._20, m._30}

	//r[0] = m[0] * n._00 + m[1] * n._10 + m[2] * n._20 + m[3] * n._30; // n[0]; n._10 = n[0][3-1] = n[0][2] => shuf<2>
	//r[1] = m[0] * n._01 + m[1] * n._11 + m[2] * n._21 + m[3] * n._31; // n[1]
	//r[2] = m[0] * n._02 + m[1] * n._12 + m[2] * n._22 + m[3] * n._32; // n[2]
	//r[3] = m[0] * n._03 + m[1] * n._13 + m[2] * n._23 + m[3] * n._33; // n[3]

	//IACA_FROM
	// 18.89, 16

	r[0] =
		_mm_add_ps(
			_mm_add_ps(
				_mm_mul_ps(m[0], _mm_shuffle_ps(n[0], n[0], _MM_SHUFFLE(3,3,3,3))),
				_mm_mul_ps(m[1], _mm_shuffle_ps(n[0], n[0], _MM_SHUFFLE(2,2,2,2)))
			),
			_mm_add_ps(
				_mm_mul_ps(m[2], _mm_shuffle_ps(n[0], n[0], _MM_SHUFFLE(1,1,1,1))),
				_mm_mul_ps(m[3], _mm_shuffle_ps(n[0], n[0], _MM_SHUFFLE(0,0,0,0)))
			)
		);

	r[1] =
		_mm_add_ps(
			_mm_add_ps(
				_mm_mul_ps(m[0], _mm_shuffle_ps(n[1], n[1], _MM_SHUFFLE(3,3,3,3))),
				_mm_mul_ps(m[1], _mm_shuffle_ps(n[1], n[1], _MM_SHUFFLE(2,2,2,2)))
			),
			_mm_add_ps(
				_mm_mul_ps(m[2], _mm_shuffle_ps(n[1], n[1], _MM_SHUFFLE(1,1,1,1))),
				_mm_mul_ps(m[3], _mm_shuffle_ps(n[1], n[1], _MM_SHUFFLE(0,0,0,0)))
			)
		);

	r[2] =
		_mm_add_ps(
			_mm_add_ps(
				_mm_mul_ps(m[0], _mm_shuffle_ps(n[2], n[2], _MM_SHUFFLE(3,3,3,3))),
				_mm_mul_ps(m[1], _mm_shuffle_ps(n[2], n[2], _MM_SHUFFLE(2,2,2,2)))
			),
			_mm_add_ps(
				_mm_mul_ps(m[2], _mm_shuffle_ps(n[2], n[2], _MM_SHUFFLE(1,1,1,1))),
				_mm_mul_ps(m[3], _mm_shuffle_ps(n[2], n[2], _MM_SHUFFLE(0,0,0,0)))
			)
		);

	r[3] =
		_mm_add_ps(
			_mm_add_ps(
				_mm_mul_ps(m[0], _mm_shuffle_ps(n[3], n[3], _MM_SHUFFLE(3,3,3,3))),
				_mm_mul_ps(m[1], _mm_shuffle_ps(n[3], n[3], _MM_SHUFFLE(2,2,2,2)))
			),
			_mm_add_ps(
				_mm_mul_ps(m[2], _mm_shuffle_ps(n[3], n[3], _MM_SHUFFLE(1,1,1,1))),
				_mm_mul_ps(m[3], _mm_shuffle_ps(n[3], n[3], _MM_SHUFFLE(0,0,0,0)))
			)
		);

	//IACA_TO
}

struct alignas(sizeof(__m128)) m128 {
	__m128 v;
	m128() {} // = default makes overhead?
	m128(__m128 const a) : v(a) {}
	operator __m128() const { return v; }
};

m128 operator + (m128 const a, m128 const b) { return _mm_add_ps(a, b); }
m128 operator - (m128 const a, m128 const b) { return _mm_sub_ps(a, b); }
m128 operator * (m128 const a, m128 const b) { return _mm_mul_ps(a, b); }
m128 operator / (m128 const a, m128 const b) { return _mm_div_ps(a, b); }

template <int a, int b, int c, int d> __m128 shuf(__m128 const u, __m128 const v)
{ return _mm_shuffle_ps(u, v, _MM_SHUFFLE(a, b, c, d)); }
template <int i> __m128 shuf(__m128 const u, __m128 const v)
{ return _mm_shuffle_ps(u, v, _MM_SHUFFLE(i, i, i, i)); }
template <int a, int b, int c, int d> __m128 shuf(__m128 const v)
{ return _mm_shuffle_ps(v, v, _MM_SHUFFLE(a, b, c, d)); }
template <int i> __m128 shuf(__m128 const v)
{ return _mm_shuffle_ps(v, v, _MM_SHUFFLE(i, i, i, i)); }
template <int a, int b, int c, int d> __m128 shufd(__m128 const v)
{ return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(v), _MM_SHUFFLE(a, b, c, d))); }
template <int i> __m128 shufd(__m128 const v)
{ return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(v), _MM_SHUFFLE(i, i, i, i))); }

void __vmathcall mul_mtx4_mtx4_sse_v2(__m128* const r, __m128 const* const m, __m128 const* const n) {
	//IACA_FROM
	// 19; 16
	r[0] = m[0]*shuf<3>(n[0]) + m[1]*shuf<2>(n[0]) + m[2]*shuf<1>(n[0]) + m[3]*shuf<0>(n[0]);
	r[1] = m[0]*shuf<3>(n[1]) + m[1]*shuf<2>(n[1]) + m[2]*shuf<1>(n[1]) + m[3]*shuf<0>(n[1]);
	r[2] = m[0]*shuf<3>(n[2]) + m[1]*shuf<2>(n[2]) + m[2]*shuf<1>(n[2]) + m[3]*shuf<0>(n[2]);
	r[3] = m[0]*shuf<3>(n[3]) + m[1]*shuf<2>(n[3]) + m[2]*shuf<1>(n[3]) + m[3]*shuf<0>(n[3]);
	//IACA_TO
}

__m128 mad(__m128 const a, __m128 const b, __m128 const c) {
	return _mm_add_ps(_mm_mul_ps(a, b), c);
}

void __vmathcall mul_mtx4_mtx4_sse_v3(__m128* const r, __m128 const* const m, __m128 const* const n) {
	//IACA_FROM
	// 18.89; 16
	r[0] = mad(m[0], shuf<3>(n[0]), m[1]*shuf<2>(n[0])) + mad(m[2], shuf<1>(n[0]), m[3]*shuf<0>(n[0]));
	r[1] = mad(m[0], shuf<3>(n[1]), m[1]*shuf<2>(n[1])) + mad(m[2], shuf<1>(n[1]), m[3]*shuf<0>(n[1]));
	r[2] = mad(m[0], shuf<3>(n[2]), m[1]*shuf<2>(n[2])) + mad(m[2], shuf<1>(n[2]), m[3]*shuf<0>(n[2]));
	r[3] = mad(m[0], shuf<3>(n[3]), m[1]*shuf<2>(n[3])) + mad(m[2], shuf<1>(n[3]), m[3]*shuf<0>(n[3]));
	//IACA_TO
}

void __vmathcall mul_mtx4_mtx4_sse_v4(__m128* const r, __m128 const* const m, __m128 const* const n) {
	//IACA_FROM
	// 18.89; 16
	_mm_stream_ps(reinterpret_cast<float*>(r + 0),
		mad(m[0], shuf<3>(n[0]), m[1]*shuf<2>(n[0])) + mad(m[2], shuf<1>(n[0]), m[3]*shuf<0>(n[0])));

	_mm_stream_ps(reinterpret_cast<float*>(r + 1),
		mad(m[0], shuf<3>(n[1]), m[1]*shuf<2>(n[1])) + mad(m[2], shuf<1>(n[1]), m[3]*shuf<0>(n[1])));

	_mm_stream_ps(reinterpret_cast<float*>(r + 2),
		mad(m[0], shuf<3>(n[2]), m[1]*shuf<2>(n[2])) + mad(m[2], shuf<1>(n[2]), m[3]*shuf<0>(n[2])));

	_mm_stream_ps(reinterpret_cast<float*>(r + 3),
		mad(m[0], shuf<3>(n[3]), m[1]*shuf<2>(n[3])) + mad(m[2], shuf<1>(n[3]), m[3]*shuf<0>(n[3])));
	//IACA_TO
}
/*
void __vmathcall mul_mtx4_mtx4_sse_v5(__m128* const _r, __m128 const* const m, __m128 const* const n) {
	mtx4& r = **reinterpret_cast<mtx4* const*>(&_r);
	IACA_FROM
	r._00 = _mm_cvtss_f32(_mm_dp_ps(m[0], n[0], 0x1f));
	r._01 = _mm_cvtss_f32(_mm_dp_ps(m[0], n[1], 0x1f));
	r._02 = _mm_cvtss_f32(_mm_dp_ps(m[0], n[2], 0x1f));
	r._03 = _mm_cvtss_f32(_mm_dp_ps(m[0], n[3], 0x1f));

	r._10 = _mm_cvtss_f32(_mm_dp_ps(m[1], n[0], 0x1f));
	r._11 = _mm_cvtss_f32(_mm_dp_ps(m[1], n[1], 0x1f));
	r._12 = _mm_cvtss_f32(_mm_dp_ps(m[1], n[2], 0x1f));
	r._13 = _mm_cvtss_f32(_mm_dp_ps(m[1], n[3], 0x1f));

	r._20 = _mm_cvtss_f32(_mm_dp_ps(m[2], n[0], 0x1f));
	r._21 = _mm_cvtss_f32(_mm_dp_ps(m[2], n[1], 0x1f));
	r._22 = _mm_cvtss_f32(_mm_dp_ps(m[2], n[2], 0x1f));
	r._23 = _mm_cvtss_f32(_mm_dp_ps(m[2], n[3], 0x1f));

	r._30 = _mm_cvtss_f32(_mm_dp_ps(m[3], n[0], 0x1f));
	r._31 = _mm_cvtss_f32(_mm_dp_ps(m[3], n[1], 0x1f));
	r._32 = _mm_cvtss_f32(_mm_dp_ps(m[3], n[2], 0x1f));
	r._33 = _mm_cvtss_f32(_mm_dp_ps(m[3], n[3], 0x1f));
	IACA_TO
}
//*/
void __vmathcall mul_mtx4_mtx4_avx_v1(__m128* const r, __m128 const* const m, __m128 const* const n) {
	/*
	r00 = m00*n00 + m01*n10 + m02*n20 + m03*n30 // m0:m1 * [n00]:[n10] + m2:m3 * [n20]:[n30]
	r10 = m10*n00 + m11*n10 + m12*n20 + m13*n30
	r20 = m20*n00 + m21*n10 + m22*n20 + m23*n30
	r30 = m30*n00 + m31*n10 + m32*n20 + m33*n30

	r01 = m00*n01 + m01*n11 + m02*n21 + m03*n31 // m0:m1 * [n01]:[n11] + m2:m3 * [n21]:[n31]
	r11 = m10*n01 + m11*n11 + m12*n21 + m13*n31
	r21 = m20*n01 + m21*n11 + m22*n21 + m23*n31
	r31 = m30*n01 + m31*n11 + m32*n21 + m33*n31

	r02 = m00*n02 + m01*n12 + m02*n22 + m03*n32 // m0:m1 * [n02]:[n12] + m2:m3 * [n22]:[n32]
	r12 = m10*n02 + m11*n12 + m12*n22 + m13*n32
	r22 = m20*n02 + m21*n12 + m22*n22 + m23*n32
	r32 = m30*n02 + m31*n12 + m32*n22 + m33*n32

	r03 = m00*n03 + m01*n13 + m02*n23 + m03*n33 // m0:m1 * [n03]:[n13] + m2:m3 * [n23]:[n33]
	r13 = m10*n03 + m11*n13 + m12*n23 + m13*n33
	r23 = m20*n03 + m21*n13 + m22*n23 + m23*n33
	r33 = m30*n03 + m31*n13 + m32*n23 + m33*n33
-
	r0 = m0*n00 + m1*n10 + m2*n20 + m3*n30
	r1 = m0*n01 + m1*n11 + m2*n21 + m3*n31
	r2 = m0*n02 + m1*n12 + m2*n22 + m3*n32
	r3 = m0*n03 + m1*n13 + m2*n23 + m3*n33
	*/

	/*
	_00, _10, _20, _30,
	_01, _11, _21, _31,
	_02, _12, _22, _32,
	_03, _13, _23, _33

	n0n1 = _00, _10, _20, _30 : _01, _11, _21, _31
	n2n3 = _02, _12, _22, _32 : _03, _13, _23, _33

	mm[0] = m0:m0
	mm[1] = m1:m1 
	mm[2] = m2:m2
	mm[3] = m3:m3

	r0r1 =
	mm[0] * n00,n00,n00,n00:n01,n01,n01,n01 + // perm<3>
	mm[1] * n10,n10,n10,n10:n11,n11,n11,n11 + // perm<2>
	mm[2] * n20,n20,n20,n20:n21,n21,n21,n21 + // perm<1>
	mm[3] * n30,n30,n30,n30:n31,n31,n31,n31   // perm<0>

	r2r3 =
	mm[0] * n02,n02,n02,n02:n03,n03,n03,n03 +
	mm[1] * n12,n12,n12,n12:n13,n13,n13,n13 +
	mm[2] * n22,n22,n22,n22:n23,n23,n23,n23 +
	mm[3] * n32,n32,n32,n32:n33,n33,n33,n33

	r0r1 = mm[0]*n0n1<3,3,3,3> + mm[1]*n0n1<2,2,2,2> + mm[2]*n0n1<1,1,1,1> + mm[3]*n0n1<0,0,0,0>
	r2r3 = mm[0]*n2n3<3,3,3,3> + mm[1]*n2n3<2,2,2,2> + mm[2]*n2n3<1,1,1,1> + mm[3]*n2n3<0,0,0,0>

	r0r1 = mm[0]*n0n1<3> + mm[1]*n0n1<2> + mm[2]*n0n1<1> + mm[3]*n0n1<0>
	r2r3 = mm[0]*n2n3<3> + mm[1]*n2n3<2> + mm[2]*n2n3<1> + mm[3]*n2n3<0>
	*/

	//IACA_FROM
	// 13; 12
	__m256 mm0 = _mm256_set_m128(m[0], m[0]);
	__m256 mm1 = _mm256_set_m128(m[1], m[1]);
	__m256 mm2 = _mm256_set_m128(m[2], m[2]);
	__m256 mm3 = _mm256_set_m128(m[3], m[3]);

	__m256 n0n1 = _mm256_load_ps(reinterpret_cast<float const*>(n + 0));
	__m256 y1 = _mm256_permute_ps(n0n1, 0xFF);//3,3,3,3; 0b11'11'11'11; 0b1111'1111
	__m256 y2 = _mm256_permute_ps(n0n1, 0xAA);//2,2,2,2; 0b10'10'10'10; 0b1010'1010
	__m256 y3 = _mm256_permute_ps(n0n1, 0x55);//1,1,1,1; 0b01'01'01'01; 0b0101'0101
	__m256 y4 = _mm256_permute_ps(n0n1, 0x00);//0,0,0,0; 0b00'00'00'00; 0b0000'0000

	y1 = _mm256_mul_ps(y1, mm0);
	y2 = _mm256_mul_ps(y2, mm1);
	y3 = _mm256_mul_ps(y3, mm2);
	y4 = _mm256_mul_ps(y4, mm3);

	y1 = _mm256_add_ps(y1, y2);
	y3 = _mm256_add_ps(y3, y4);
	y1 = _mm256_add_ps(y1, y3);

	__m256 n2n3 = _mm256_load_ps(reinterpret_cast<float const*>(n + 2));
	__m256 y5 = _mm256_permute_ps(n2n3, 0xFF);
	__m256 y6 = _mm256_permute_ps(n2n3, 0xAA);
	__m256 y7 = _mm256_permute_ps(n2n3, 0x55);
	__m256 y8 = _mm256_permute_ps(n2n3, 0x00);

	y5 = _mm256_mul_ps(y5, mm0);
	y6 = _mm256_mul_ps(y6, mm1);
	y7 = _mm256_mul_ps(y7, mm2);
	y8 = _mm256_mul_ps(y8, mm3);

	y5 = _mm256_add_ps(y5, y6);
	y7 = _mm256_add_ps(y7, y8);
	y5 = _mm256_add_ps(y5, y7);

	_mm256_stream_ps(reinterpret_cast<float* const>(r)+0, y1);
	_mm256_stream_ps(reinterpret_cast<float* const>(r)+8, y5);
	//IACA_TO
}

void __vmathcall mul_mtx4_mtx4_avx_v1m(__m128* const _r, __m128 const* const m, __m128 const* const n) {
	mtx4& r = *reinterpret_cast<mtx4* const>(_r);

	//IACA_FROM

	__m256 mm0 = _mm256_set_m128(m[0], m[0]);
	__m256 mm1 = _mm256_set_m128(m[1], m[1]);
	__m256 mm2 = _mm256_set_m128(m[2], m[2]);
	__m256 mm3 = _mm256_set_m128(m[3], m[3]);

	__m256 n0n1 = _mm256_load_ps(reinterpret_cast<float const*>(n + 0));
	__m256 y1 = _mm256_permute_ps(n0n1, 0xFF);//3,3,3,3; 0b11'11'11'11; 0b1111'1111
	__m256 y2 = _mm256_permute_ps(n0n1, 0xAA);//2,2,2,2; 0b10'10'10'10; 0b1010'1010
	__m256 y3 = _mm256_permute_ps(n0n1, 0x55);//1,1,1,1; 0b01'01'01'01; 0b0101'0101
	__m256 y4 = _mm256_permute_ps(n0n1, 0x00);//0,0,0,0; 0b00'00'00'00; 0b0000'0000

	y1 = _mm256_mul_ps(y1, mm0);
	y2 = _mm256_mul_ps(y2, mm1);
	y3 = _mm256_mul_ps(y3, mm2);
	y4 = _mm256_mul_ps(y4, mm3);

	y1 = _mm256_add_ps(y1, y2);
	y3 = _mm256_add_ps(y3, y4);
	r.ymm[0] = _mm256_add_ps(y1, y3);

	__m256 n2n3 = _mm256_load_ps(reinterpret_cast<float const*>(n + 2));
	__m256 y5 = _mm256_permute_ps(n2n3, 0xFF);
	__m256 y6 = _mm256_permute_ps(n2n3, 0xAA);
	__m256 y7 = _mm256_permute_ps(n2n3, 0x55);
	__m256 y8 = _mm256_permute_ps(n2n3, 0x00);

	y5 = _mm256_mul_ps(y5, mm0);
	y6 = _mm256_mul_ps(y6, mm1);
	y7 = _mm256_mul_ps(y7, mm2);
	y8 = _mm256_mul_ps(y8, mm3);

	y5 = _mm256_add_ps(y5, y6);
	y7 = _mm256_add_ps(y7, y8);
	r.ymm[1] = _mm256_add_ps(y5, y7);

	//IACA_TO
}

struct alignas(sizeof(__m256)) m256 {
	__m256 v;
	m256() {}
	m256(__m256 const a) : v(a) {}
	operator __m256() const { return v; }
};

m256 operator + (m256 const a, m256 const b) { return _mm256_add_ps(a, b); }
m256 operator - (m256 const a, m256 const b) { return _mm256_sub_ps(a, b); }
m256 operator * (m256 const a, m256 const b) { return _mm256_mul_ps(a, b); }
m256 operator / (m256 const a, m256 const b) { return _mm256_div_ps(a, b); }

template <int i> __m256 perm(__m256 const v)
{ return _mm256_permute_ps(v, _MM_SHUFFLE(i, i, i, i)); }
template <int a, int b, int c, int d> __m256 perm(__m256 const v)
{ return _mm256_permute_ps(v, _MM_SHUFFLE(a, b, c, d)); }
template <int i, int j> __m256 perm(__m256 const v)
{ return _mm256_permutevar_ps(v, _mm256_set_epi32(i, i, i, i, j, j, j, j)); }
template <int a, int b, int c, int d, int e, int f, int g, int h> __m256 perm(__m256 const v)
{ return _mm256_permutevar_ps(v, _mm256_set_epi32(a, b, c, d, e, f, g, h)); }

__m256 mad(__m256 const a, __m256 const b, __m256 const c) {
	return _mm256_add_ps(_mm256_mul_ps(a, b), c);
}

void __vmathcall mul_mtx4_mtx4_avx_v2(__m128* const r, __m128 const* const m, __m128 const* const n) {
	//IACA_FROM
	// 10; 8.58
	//broadcast makes some profit
	__m256 const mm[] {
		_mm256_broadcast_ps(m+0),
		_mm256_broadcast_ps(m+1),
		_mm256_broadcast_ps(m+2),
		_mm256_broadcast_ps(m+3)
	};

	__m256 const n0n1 = _mm256_load_ps(reinterpret_cast<float const*>(n + 0));
	_mm256_stream_ps(reinterpret_cast<float*>(r),
		mad(perm<3>(n0n1), mm[0], perm<2>(n0n1)*mm[1])+mad(perm<1>(n0n1), mm[2], perm<0>(n0n1)*mm[3]));

	__m256 const n2n3 = _mm256_load_ps(reinterpret_cast<float const*>(n + 2));
	_mm256_stream_ps(reinterpret_cast<float*>(r)+8,
		mad(perm<3>(n2n3), mm[0], perm<2>(n2n3)*mm[1])+mad(perm<1>(n2n3), mm[2], perm<0>(n2n3)*mm[3]));

	//IACA_TO
}

void __vmathcall mul_mtx4_mtx4_avx_v2m(__m128* const _r, __m128 const* const m, __m128 const* const n) {
	mtx4& r = *reinterpret_cast<mtx4* const>(_r);
	
	//IACA_FROM
	
	__m256 const mm[] {
		_mm256_broadcast_ps(m+0),
		_mm256_broadcast_ps(m+1),
		_mm256_broadcast_ps(m+2),
		_mm256_broadcast_ps(m+3)
	};

	__m256 const n0n1 = _mm256_load_ps(reinterpret_cast<float const*>(n + 0));
	r.ymm[0] = mad(perm<3>(n0n1), mm[0], perm<2>(n0n1)*mm[1])+mad(perm<1>(n0n1), mm[2], perm<0>(n0n1)*mm[3]);

	__m256 const n2n3 = _mm256_load_ps(reinterpret_cast<float const*>(n + 2));
	r.ymm[1] = mad(perm<3>(n2n3), mm[0], perm<2>(n2n3)*mm[1])+mad(perm<1>(n2n3), mm[2], perm<0>(n2n3)*mm[3]);

	//IACA_TO
}

__m256 fma(__m256 const a, __m256 const b, __m256 const c) {
	return _mm256_fmadd_ps(a, b, c);
}

void __vmathcall mul_mtx4_mtx4_avx_fma(__m128* const r, __m128 const* const m, __m128 const* const n) {
	//IACA_FROM
	// 9.21; 8

	__m256 const mm[]{
		_mm256_broadcast_ps(m + 0),
		_mm256_broadcast_ps(m + 1),
		_mm256_broadcast_ps(m + 2),
		_mm256_broadcast_ps(m + 3) };

	__m256 const n0n1 = _mm256_load_ps(reinterpret_cast<float const*>(n + 0));
	_mm256_stream_ps(reinterpret_cast<float*>(r + 0),
		fma(perm<3>(n0n1), mm[0], perm<2>(n0n1)*mm[1])+fma(perm<1>(n0n1), mm[2], perm<0>(n0n1)*mm[3]));

	__m256 const n2n3 = _mm256_load_ps(reinterpret_cast<float const*>(n + 2));
	_mm256_stream_ps(reinterpret_cast<float*>(r + 2),
		fma(perm<3>(n2n3), mm[0], perm<2>(n2n3)*mm[1])+fma(perm<1>(n2n3), mm[2], perm<0>(n2n3)*mm[3]));

	//IACA_TO
}

void __vmathcall mul_mtx4_mtx4_avx_fma_m(__m128* const _r, __m128 const* const m, __m128 const* const n) {
	mtx4& r = **reinterpret_cast<mtx4* const*>(&_r);

	//###
	IACA_FROM
	__m256 const mm[]{
		_mm256_broadcast_ps(m + 0),
		_mm256_broadcast_ps(m + 1),
		_mm256_broadcast_ps(m + 2),
		_mm256_broadcast_ps(m + 3) };

	__m256 const n0n1 = _mm256_load_ps(reinterpret_cast<float const*>(n + 0));
	r.ymm[0] = fma(perm<3>(n0n1), mm[0], perm<2>(n0n1)*mm[1])+fma(perm<1>(n0n1), mm[2], perm<0>(n0n1)*mm[3]);

	__m256 const n2n3 = _mm256_load_ps(reinterpret_cast<float const*>(n + 2));
	r.ymm[1] = fma(perm<3>(n2n3), mm[0], perm<2>(n2n3)*mm[1])+fma(perm<1>(n2n3), mm[2], perm<0>(n2n3)*mm[3]);
	IACA_TO
}

#ifdef __AVX512F__

struct alignas(sizeof(__m512)) m512 {
	__m512 v;
	m512() {}
	m512(__m512 const a) : v(a) {}
	operator __m512() const { return v; }
};

m512 operator + (m512 const a, m512 const b) { return _mm512_add_ps(a, b); }
m512 operator - (m512 const a, m512 const b) { return _mm512_sub_ps(a, b); }
m512 operator * (m512 const a, m512 const b) { return _mm512_mul_ps(a, b); }
m512 operator / (m512 const a, m512 const b) { return _mm512_div_ps(a, b); }

template <int i> __m512 perm(__m512 const v)
{ return _mm512_permute_ps(v, _MM_SHUFFLE(i, i, i, i)); }
template <int a, int b, int c, int d> __m512 perm(__m512 const v)
{ return _mm512_permute_ps(v, _MM_SHUFFLE(a, b, c, d)); }

__m512 fma(__m512 const a, __m512 const b, __m512 const c) {
	return _mm512_fmadd_ps(a, b, c);
}

void __vmathcall mul_mtx4_mtx4_avx512(__m128* const r, __m128 const* const m, __m128 const* const _n) {
	//IACA_FROM
	// 4.79; 5.42
	__m512 const mm[]{
		_mm512_broadcast_f32x4(m[0]),
		_mm512_broadcast_f32x4(m[1]),
		_mm512_broadcast_f32x4(m[2]),
		_mm512_broadcast_f32x4(m[3]) };

	__m512 const n = _mm512_load_ps(reinterpret_cast<float const*>(_n));
	_mm512_stream_ps(reinterpret_cast<float*>(r),
		fma(perm<3>(n), mm[0], perm<2>(n)*mm[1])+fma(perm<1>(n), mm[2], perm<0>(n)*mm[3]));
	//IACA_TO
}

void __vmathcall mul_mtx4_mtx4_avx512m(__m128* const _r, __m128 const* const m, __m128 const* const _n) {
	mtx4& r = **reinterpret_cast<mtx4* const*>(&_r);
	//IACA_FROM
	__m512 const mm[]{
		_mm512_broadcast_f32x4(m[0]),
		_mm512_broadcast_f32x4(m[1]),
		_mm512_broadcast_f32x4(m[2]),
		_mm512_broadcast_f32x4(m[3]) };

	__m512 const n = _mm512_load_ps(reinterpret_cast<float const*>(_n + 0));
	r.zmm = fma(perm<3>(n), mm[0], perm<2>(n)*mm[1])+fma(perm<1>(n), mm[2], perm<0>(n)*mm[3]);
	//IACA_TO
}

#endif // __AVX512F__

mtx4 get_mtx4_1() {
	return mtx4(
		1.f, 2.f, 3.f, 4.f,
		5.f, 6.f, 7.f, 8.f,
		4.f, 3.f, 2.f, 1.f,
		8.f, 7.f, 6.f, 5.f);
}

mtx4 get_mtx4_2() {
	return mtx4(
		3.f, 7.f, 2.f, 4.f,
		6.f, 3.f, 5.f, 1.f,
		9.f, 6.f, 8.f, 2.f,
		4.f, 2.f, 7.f, 6.f);
}

using ts = rts;

template <typename F>
ts speed_test(F&& f, mtx4 const* const src, mtx4* const dst, std::size_t const size, std::size_t const icount) {
	mtx4 const identity = mtx4::identity();

	ts t;
	for (size_t c = 0; c < icount; c++) {
		for (size_t i = 0; i < size; i++) {
			f(dst[i], src[i], identity);
		}
	}
	t.stop();

	return t;
}

void mul_mtx4_mtx4_idle(__m128* const, __m128 const* const, __m128 const* const) {
#	ifdef WIN32
	__nop();
#	endif
}

#define SPEED_TEST(f, idx) \
	{ \
		mtx4 const identity = mtx4::identity(); \
		ts t; \
		for (size_t c = 0; c < icount; ++c) { \
			for (size_t i = 0; i < size; ++i) { \
				f(dst[i], src[i], identity); \
			} \
		} \
		t.stop(); \
		t -= t_idle.passed(); \
		tss[idx] = t; \
	}

int main() {
	static_assert(sizeof(vec4)==sizeof(float)* 4, "vec4 size is not 16");
	static_assert(sizeof(mtx4)==sizeof(float)*16, "mtx4 size is not 64");

	mtx4 mref = mtx4::zero();
	mul_mtx4_mtx4_unroll(mref, mtx4::identity(), mtx4::identity());
	assert(mref == mtx4::identity());

	mtx4 r = mtx4::zero(), m = get_mtx4_1(), n = get_mtx4_2();
	mtx4 const identity = mtx4::identity();

	mul_mtx4_mtx4_unroll(mref, m, n);

	size_t const icount = 10'000;
	size_t const size   = 1000;
	mtx4* src = (mtx4*) _mm_malloc(sizeof(mtx4)*size, 64);
	mtx4* dst = (mtx4*) _mm_malloc(sizeof(mtx4)*size, 64);

	for (size_t i = 0; i < size; ++i) {
		src[i] = mtx4::identity();
		dst[i] = mtx4::identity();
	}

	ts t_idle;
	for (size_t c = 0; c < icount; ++c) {
		for (size_t i = 0; i < size; ++i) {
			mul_mtx4_mtx4_idle(dst[i], src[i], identity);
		}
	}
	t_idle.stop();

	std::cout << "cpuid+cycle tics: " << t_idle.passed() << std::endl;
	std::cout << "cpuid+cycle tics per iteration: " << t_idle.passed() / (1. * size * icount) << std::endl;

	//
	//SetThreadAffinityMask(GetCurrentThread(), 1 << (0));
	//SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);

	//

	using mul_mtx4_mtx4_t = void (__vmathcall *)(__m128* const r, __m128 const* const m, __m128 const* const n);
	using func_data_t = std::tuple<std::string, mul_mtx4_mtx4_t, float, float>;

	func_data_t funcs[] {
		{ "unroll" , mul_mtx4_mtx4_unroll    ,  69.95f,  64.00f },
		{ "loop"   , mul_mtx4_mtx4_loop      , 233.60f, 190.08f },
#       if defined GLM_VERSION
		{ "glm"    , mul_mtx4_mtx4_glm        ,  69.95f,  64.00f }, // expecting same as unroll timings
		{ "glm_sse", mul_mtx4_mtx4_glm_sse    ,  18.89f,  16.00f }, // expecting same as raw sse timings
#       endif
		{ "sse_v1" , mul_mtx4_mtx4_sse_v1    ,  18.89f,  16.00f },
		{ "sse_v2" , mul_mtx4_mtx4_sse_v2    ,  19.00f,  16.00f },
		{ "sse_v3" , mul_mtx4_mtx4_sse_v3    ,  18.89f,  16.00f },
		{ "sse_v4s", mul_mtx4_mtx4_sse_v4    ,  18.89f,  16.00f },
#		ifdef __AVX__
		{ "avx_v1m", mul_mtx4_mtx4_avx_v1m   ,  13.00f,  12.00f },
		{ "avx_v1s", mul_mtx4_mtx4_avx_v1    ,  13.00f,  12.00f },
		{ "avx_v2m", mul_mtx4_mtx4_avx_v2m   ,  10.00f,   8.58f },
		{ "avx_v2s", mul_mtx4_mtx4_avx_v2    ,  10.00f,   8.58f },
#		endif
#		ifdef __FMA__
		{ "AVX+FMAm", mul_mtx4_mtx4_avx_fma_m,   9.21f,   8.00f },
		{ "AVX+FMAs", mul_mtx4_mtx4_avx_fma  ,   9.21f,   8.00f },
#		endif
#		ifdef __AVX512F__
		{ "AVX512m", mul_mtx4_mtx4_avx512m   ,   4.79f,   5.42f  },
		{ "AVX512s", mul_mtx4_mtx4_avx512    ,   4.79f,   5.42f  },
#		endif
	};
	
	{
		ts t;
		for (size_t c = 0; c < icount; ++c) {
			for (size_t i = 0; i < size; ++i) {
				mul_mtx4_mtx4_avx_v2m(dst[i], src[i], identity);
			}
		}
		t.stop();

		std::cout << " direct call tics: " << t.passed() * (1.f / (size*icount)) << std::endl;
	}

	{
		mul_mtx4_mtx4_t f = mul_mtx4_mtx4_avx_v2m;
		ts t;
		for (size_t c = 0; c < icount; ++c) {
			for (size_t i = 0; i < size; ++i) {
				f(dst[i], src[i], identity);
			}
		}
		t.stop();

		std::cout << "pointer call tics: "  << t.passed() * (1.f / (size*icount)) << std::endl;
	}


	// check result correctness
	for (size_t i = 0; i < std::size(funcs); i++) {
#if defined GLM_VERSION
		if (std::get<1>(funcs[i]) == mul_mtx4_mtx4_glm) {
			// make sure glm::mat4 format was understood correctly
			glm::mat4 _m = *reinterpret_cast<glm::mat4 const* const>(&m);
			glm::mat4 _n = *reinterpret_cast<glm::mat4 const* const>(&n);
			for (int i = 0; i < 4; _m[i] = _m[i].wzyx(), i++);
			for (int i = 0; i < 4; _n[i] = _n[i].wzyx(), i++);
			glm::mat4 _r = _m * _n;
			for (int i = 0; i < 4; _r[i] = _r[i].wzyx(), i++);
			r = *reinterpret_cast<mtx4 const* const>(&_r);
		} else if (std::get<1>(funcs[i]) == mul_mtx4_mtx4_glm_sse) {
			// make sure glm::mat4 format was understood correctly
			glm_vec4 _m[4]{ shuf<0,1,2,3>(m[0]), shuf<0,1,2,3>(m[1]), shuf<0,1,2,3>(m[2]), shuf<0,1,2,3>(m[3]) };
			glm_vec4 _n[4]{ shuf<0,1,2,3>(n[0]), shuf<0,1,2,3>(n[1]), shuf<0,1,2,3>(n[2]), shuf<0,1,2,3>(n[3]) };
			glm_vec4 _r[4];
			mul_mtx4_mtx4_glm_sse(_r, _m, _n);
			for (int i = 0; i < 4; r[i] = shuf<0, 1, 2, 3>(_r[i]), i++);
		} else {
#endif
			std::get<1>(funcs[i])(r, m, n);
#if defined GLM_VERSION
		}
#endif
		if (!(r == mref)) {
			std::cout << i << ": " << std::get<0>(funcs[i]) << std::endl;
			std::cout << r << std::endl;
			std::cout << mref << std::endl;
		}
		assert(r == mref);
	}

	std::array<ts, std::size(funcs)> tss;

	/* slower, +70-80% tics per `mul_mtx4_mtx4_t` function call
	for (size_t i = 0; i < std::size(funcs); i++) {
		tss[i] = speed_test(std::get<1>(funcs[i]), src, dst, size, icount);
		tss[i] -= t_idle.passed();
	}
	//*/

	size_t counter = 0;
	// run speed tests
	SPEED_TEST(mul_mtx4_mtx4_unroll   , counter++)
	SPEED_TEST(mul_mtx4_mtx4_loop     , counter++);
#   if defined GLM_VERSION
	SPEED_TEST(mul_mtx4_mtx4_glm      , counter++);
	SPEED_TEST(mul_mtx4_mtx4_glm_sse  , counter++);
#   endif
	SPEED_TEST(mul_mtx4_mtx4_sse_v1   , counter++);
	SPEED_TEST(mul_mtx4_mtx4_sse_v2   , counter++);
	SPEED_TEST(mul_mtx4_mtx4_sse_v3   , counter++);
	SPEED_TEST(mul_mtx4_mtx4_sse_v4   , counter++);
#	ifdef __AVX__
	SPEED_TEST(mul_mtx4_mtx4_avx_v1m  , counter++);
	SPEED_TEST(mul_mtx4_mtx4_avx_v1   , counter++);
	SPEED_TEST(mul_mtx4_mtx4_avx_v2m  , counter++);
	SPEED_TEST(mul_mtx4_mtx4_avx_v2   , counter++);
#	endif
#	ifdef __FMA__
	SPEED_TEST(mul_mtx4_mtx4_avx_fma_m, counter++);
	SPEED_TEST(mul_mtx4_mtx4_avx_fma  , counter++);
#	endif
#	ifdef __AVX512F__
	SPEED_TEST(mul_mtx4_mtx4_avx512m  , counter++);
	SPEED_TEST(mul_mtx4_mtx4_avx512   , counter++);
#	endif

	// display results
	std::cout << "-----------------------------------------" << std::endl;

	std::cout << "      name: rcoef   coef   rtics    tics" << std::endl;

	static const double K = 1. / (size * icount);
	for (size_t i = 0; i < std::size(funcs); i++) {
		std::cout << std::setprecision(2) << std::fixed
			<< std::setw(10) << std::get<0>(funcs[i]).c_str()
			<< ": " << std::setw(5) << std::get<2>(funcs[0])/std::get<2>(funcs[i])
			<< "; " << std::setw(5) << static_cast<double>(tss[0].passed())/tss[i].passed()
			<< "; " << std::setw(6) << std::get<2>(funcs[i])
			<< "; " << std::setw(6) << std::fixed << tss[i].passed() * K
			<< std::endl;
	}

	std::cout << "-----------------------------------------" << std::endl;

	std::cout << "rcoef: performance coefficient from reference theoretical cpu tics" << std::endl;
	std::cout << "coef : performance coefficient from real cpu tics" << std::endl;
	std::cout << "rtics: reference theoretical cpu tics" << std::endl;
	std::cout << "tics : real cpu tics" << std::endl;

	std::cout << "reference value is taken from non-simd unrolled calculations" << std::endl;
	std::cout << "coef = value / reference_value, so coef > 1 means function is faster than the reference" << std::endl;

	std::cout << dst[0] << std::endl;

	_mm_free(src);
	_mm_free(dst);

	//_getch();
	return 0;
}

/*
i7-3770
x86:
    unroll:  1.00;  1.00;     70;  50.75
      loop:  0.30;  0.43; 233.60; 119.21
    sse_v1:  3.70;  1.84;  18.89;  27.51
    sse_v2:  3.68;  1.84;  19.00;  27.61
    sse_v3:  3.70;  1.86;  18.89;  27.22
   sse_v4s:  3.70;  1.87;  18.89;  27.18
   avx_v1m:  5.38;  2.64;  13.00;  19.21
   avx_v1s:  5.38;  2.53;  13.00;  20.03
   avx_v2m:  6.99;  3.93;  10.00;  12.91
   avx_v2s:  6.99;  2.93;  10.00;  17.34

x64:
    unroll:  1.00;  1.00;     70;  68.60
      loop:  0.30;  0.57; 233.60; 119.37
    sse_v1:  3.70;  3.12;  18.89;  21.98
    sse_v2:  3.68;  3.25;  19.00;  21.09
    sse_v3:  3.70;  3.09;  18.89;  22.19
   sse_v4s:  3.70;  3.06;  18.89;  22.39
   avx_v1m:  5.38;  7.13;  13.00;   9.61
   avx_v1s:  5.38;  4.06;  13.00;  16.90
   avx_v2m:  6.99;  7.45;  10.00;   9.20
   avx_v2s:  6.99;  4.68;  10.00;  14.64

i7-8700K:
x86:
    unroll:  1.00;  1.00;  69.95;  40.25
      loop:  0.30;  0.51; 233.60;  79.49
    sse_v1:  3.70;  2.09;  18.89;  19.31
    sse_v2:  3.68;  2.01;  19.00;  19.98
    sse_v3:  3.70;  2.04;  18.89;  19.69
   sse_v4s:  3.70;  2.05;  18.89;  19.67
   avx_v1m:  5.38;  2.83;  13.00;  14.22
   avx_v1s:  5.38;  2.85;  13.00;  14.13
   avx_v2m:  6.99;  3.43;  10.00;  11.73
   avx_v2s:  6.99;  3.41;  10.00;  11.81
  AVX+FMAm:  7.60;  3.88;   9.21;  10.38
  AVX+FMAs:  7.60;  3.90;   9.21;  10.32

x64:
    unroll:  1.00;  1.00;  69.95;  57.11
      loop:  0.30;  0.75; 233.60;  75.73
    sse_v1:  3.70;  3.61;  18.89;  15.83
    sse_v2:  3.68;  3.32;  19.00;  17.22
    sse_v3:  3.70;  3.59;  18.89;  15.92
   sse_v4s:  3.70;  3.53;  18.89;  16.18
   avx_v1m:  5.38;  8.12;  13.00;   7.03
   avx_v1s:  5.38;  4.40;  13.00;  12.98
   avx_v2m:  6.99; 10.57;  10.00;   5.40
   avx_v2s:  6.99;  5.01;  10.00;  11.39
  AVX+FMAm:  7.60;  5.87;   9.21;   9.73
  AVX+FMAs:  7.60;  5.82;   9.21;   9.81
*/
