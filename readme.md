# mul_m4_m4

Текстовый проект для кода статьи:
https://habr.com/ru/post/418247/

godbolt.org usage:
	copy paste cpp file content to godbolt.org
	clang 64: -std=c++20 -O3 -mavx2 -mfma -mavx512f -DIACA_MARKS_OFF
	icc 64  : -std=c++20 -O3 -mavx2 -mfma -mavx512f -DIACA_MARKS_OFF
	gcc 64  : -std=c++20 -O3 -mavx2 -mfma -DIACA_MARKS_OFF

msvc 2019 i7-10700:
cpuid+cycle tics: 10294668
cpuid+cycle tics per iteration: 1.02947
 direct call tics: 13.4842
pointer call tics: 14.7822
-----------------------------------------
      name: rcoef   coef   rtics    tics
    unroll:  1.00;  1.00;  69.95;  74.34
      loop:  0.30;  0.73; 233.60; 102.33
       glm:  1.00;  0.95;  69.95;  78.05
   glm sse:  3.70;  5.72;  18.89;  13.00
   sse v1 :  3.70;  5.72;  18.89;  13.00
   sse v2 :  3.68;  3.72;  19.00;  19.97
   sse v3 :  3.70;  5.72;  18.89;  13.00
   sse v4s:  3.70;  5.71;  18.89;  13.02
   avx v1m:  5.38;  5.31;  13.00;  14.00
   avx v1s:  5.38;  5.04;  13.00;  14.75
   avx v2m:  6.99;  5.68;  10.00;  13.08
   avx v2s:  6.99;  5.66;  10.00;  13.13
-----------------------------------------
rcoef: performance coefficient from reference theoretical cpu tics
coef : performance coefficient from real cpu tics
rtics: reference theoretical cpu tics
tics : real cpu tics
reference value is taken from non-simd unrolled calculations
coef = value / reference_value, so coef > 1 means function is faster than the reference
|  1.000,  0.000,  0.000,  0.000|
|  0.000,  1.000,  0.000,  0.000|
|  0.000,  0.000,  1.000,  0.000|
|  0.000,  0.000,  0.000,  1.000|


clang 64 trunk c++20 out:
cpuid+cycle tics: 482
cpuid+cycle tics per iteration: 4.82e-05
 direct call tics: 8.23507
pointer call tics: 8.00988
-----------------------------------------
      name: rcoef   coef   rtics    tics
    unroll:  1.00;  1.00;  69.95;  76.37
      loop:  0.30;  0.35; 233.60; 216.74
       glm:  1.00;  3.09;  69.95;  24.73
   glm sse:  3.70;  5.73;  18.89;  13.34
   sse v1 :  3.70;  5.78;  18.89;  13.22
   sse v2 :  3.68;  5.46;  19.00;  13.98
   sse v3 :  3.70;  5.76;  18.89;  13.25
   sse v4s:  3.70;  2.31;  18.89;  33.11
   avx v1m:  5.38; 10.27;  13.00;   7.44
   avx v1s:  5.38;  0.43;  13.00; 178.78
   avx v2m:  6.99; 10.26;  10.00;   7.45
   avx v2s:  6.99;  2.63;  10.00;  29.02
  AVX+FMAm:  7.60; 12.20;   9.21;   6.26
  AVX+FMAs:  7.60;  2.60;   9.21;  29.33
   AVX512m: 14.60; 14.12;   4.79;   5.41
   AVX512s: 14.60;  2.71;   4.79;  28.14
-----------------------------------------
rcoef: performance coefficient from reference theoretical cpu tics
coef : performance coefficient from real cpu tics
rtics: reference theoretical cpu tics
tics : real cpu tics
reference value is taken from non-simd unrolled calculations
coef = value / reference_value, so coef > 1 means function is faster than the reference
|  1.000,  0.000,  0.000,  0.000|
|  0.000,  1.000,  0.000,  0.000|
|  0.000,  0.000,  1.000,  0.000|
|  0.000,  0.000,  0.000,  1.000|

gcc 64 trunk c++20 out:
cpuid+cycle tics: 36
cpuid+cycle tics per iteration: 3.6e-06
 direct call tics: 6.47379
pointer call tics: 6.33078
-----------------------------------------
      name: rcoef   coef   rtics    tics
    unroll:  1.00;  1.00;  69.95;  33.62
      loop:  0.30;  1.07; 233.60;  31.41
       glm:  1.00;  0.80;  69.95;  41.90
   glm_sse:  3.70;  2.11;  18.89;  15.93
   sse v1 :  3.70;  0.20;  18.89; 166.38
   sse v2 :  3.68;  1.45;  19.00;  23.19
   sse v3 :  3.70;  1.30;  18.89;  25.94
   sse v4s:  3.70;  0.85;  18.89;  39.38
   avx v1m:  5.38;  1.78;  13.00;  18.85
   avx v1s:  5.38;  1.00;  13.00;  33.60
   avx v2m:  6.99;  0.20;  10.00; 164.19
   avx v2s:  6.99;  1.02;  10.00;  32.86
  AVX+FMAm:  7.60;  2.31;   9.21;  14.55
  AVX+FMAs:  7.60;  0.98;   9.21;  34.14
-----------------------------------------
rcoef: performance coefficient from reference theoretical cpu tics
coef : performance coefficient from real cpu tics
rtics: reference theoretical cpu tics
tics : real cpu tics
reference value is taken from non-simd unrolled calculations
coef = value / reference_value, so coef > 1 means function is faster than the reference
|  1.000,  0.000,  0.000,  0.000|
|  0.000,  1.000,  0.000,  0.000|
|  0.000,  0.000,  1.000,  0.000|
|  0.000,  0.000,  0.000,  1.000|

x86-64 icc 2021.2.0
cpuid+cycle tics: 346
cpuid+cycle tics per iteration: 3.46e-05
 direct call tics: 7.60582
pointer call tics: 10.2304
-----------------------------------------
      name: rcoef   coef   rtics    tics
    unroll:  1.00;  1.00;  69.95;  37.63
      loop:  0.30;  0.99; 233.60;  37.97
       glm:  1.00;  0.16;  69.95; 229.58
   glm sse:  3.70;  3.63;  18.89;  10.37
   sse v1 :  3.70;  2.50;  18.89;  15.02
   sse v2 :  3.68;  2.74;  19.00;  13.75
   sse v3 :  3.70;  2.51;  18.89;  14.99
   sse v4s:  3.70;  1.11;  18.89;  33.88
   avx v1m:  5.38;  3.19;  13.00;  11.79
   avx v1s:  5.38;  1.25;  13.00;  30.04
   avx v2m:  6.99;  0.24;  10.00; 158.60
   avx v2s:  6.99;  1.26;  10.00;  29.93
  AVX+FMAm:  7.60;  5.05;   9.21;   7.45
  AVX+FMAs:  7.60;  1.26;   9.21;  29.82
   AVX512m: 14.60;  5.59;   4.79;   6.73
   AVX512s: 14.60;  1.28;   4.79;  29.42
-----------------------------------------
rcoef: performance coefficient from reference theoretical cpu tics
coef : performance coefficient from real cpu tics
rtics: reference theoretical cpu tics
tics : real cpu tics
reference value is taken from non-simd unrolled calculations
coef = value / reference_value, so coef > 1 means function is faster than the reference
|  1.000,  0.000,  0.000,  0.000|
|  0.000,  1.000,  0.000,  0.000|
|  0.000,  0.000,  1.000,  0.000|
|  0.000,  0.000,  0.000,  1.000|

KNOWN ISSUES:
For some reason gcc discards prefix _mm512_broadcast's for AVX512 functions

gcc trunk x86-64
mul_mtx4_mtx4_avx512m(float __vector(4)*, float __vector(4) const*, float __vector(4) const*):
  vmovaps zmm0, ZMMWORD PTR [rdx]
  vpermilps zmm2, zmm0, 170
  vpermilps zmm1, zmm0, 255
  vmulps zmm2, zmm2, XMMWORD PTR [rsi+16]{1to16}
  vfmadd132ps zmm1, zmm2, XMMWORD PTR [rsi]{1to16}
  vpermilps zmm2, zmm0, 0
  vpermilps zmm0, zmm0, 85
  vmulps zmm2, zmm2, XMMWORD PTR [rsi+48]{1to16}
  vfmadd132ps zmm0, zmm2, XMMWORD PTR [rsi+32]{1to16}
  vaddps zmm0, zmm1, zmm0
  vmovaps ZMMWORD PTR [rdi], zmm0
  vzeroupper
  ret

clang
mul_mtx4_mtx4_avx512m(float __vector(4)*, float __vector(4) const*, float __vector(4) const*): # @mul_mtx4_mtx4_avx512m(float __vector(4)*, float __vector(4) const*, float __vector(4) const*)
  vbroadcastf32x4 zmm0, xmmword ptr [rsi] # zmm0 = mem[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
  vbroadcastf32x4 zmm1, xmmword ptr [rsi + 16] # zmm1 = mem[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
  vbroadcastf32x4 zmm2, xmmword ptr [rsi + 32] # zmm2 = mem[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
  vbroadcastf32x4 zmm3, xmmword ptr [rsi + 48] # zmm3 = mem[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]
  vmovaps zmm4, zmmword ptr [rdx]
  vpermilps zmm5, zmm4, 255 # zmm5 = zmm4[3,3,3,3,7,7,7,7,11,11,11,11,15,15,15,15]
  vpermilps zmm6, zmm4, 170 # zmm6 = zmm4[2,2,2,2,6,6,6,6,10,10,10,10,14,14,14,14]
  vmulps zmm1, zmm1, zmm6
  vfmadd231ps zmm1, zmm0, zmm5 # zmm1 = (zmm0 * zmm5) + zmm1
  vpermilps zmm0, zmm4, 85 # zmm0 = zmm4[1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13]
  vpermilps zmm4, zmm4, 0 # zmm4 = zmm4[0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12]
  vmulps zmm3, zmm3, zmm4
  vfmadd231ps zmm3, zmm2, zmm0 # zmm3 = (zmm2 * zmm0) + zmm3
  vaddps zmm0, zmm1, zmm3
  vmovaps zmmword ptr [rdi], zmm0
  vzeroupper
  ret
