/* -----------------------------------------------------------------------------
The copyright in this software is being made available under the BSD
License, included below. No patent rights, trademark rights and/or
other Intellectual Property Rights other than the copyrights concerning
the Software are granted under this license.

For any license concerning other Intellectual Property rights than the software,
especially patent licenses, a separate Agreement needs to be closed.
For more information please contact:

Fraunhofer Heinrich Hertz Institute
Einsteinufer 37
10587 Berlin, Germany
www.hhi.fraunhofer.de/vvc
vvc@hhi.fraunhofer.de

Copyright (c) 2018-2020, Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Fraunhofer nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.


------------------------------------------------------------------------------------------- */

/** \file     SampleAdaptiveOffsetX86.h
    \brief    SAO filter class
*/
#include "CommonDefX86.h"
#include "../SampleAdaptiveOffset.h"

//! \ingroup CommonLib
//! \{

#ifdef TARGET_SIMD_X86
#  if defined _MSC_VER
#    include <tmmintrin.h>
#  else
#    include <immintrin.h>
#  endif

#  define SAO_NUM_OFFSETS 4                           /* number of SAO offset values */
#  define SAO_EO_NUM_CATEGORIES (SAO_NUM_OFFSETS + 1) /* number of different eo categories */

#  if USE_AVX2 && !defined(_mm256_set_m128i)
#    define VVCLIB_OWN_mm256_set_m128i
#    define _mm256_set_m128i(v0, v1) _mm256_inserti128_si256(_mm256_castsi128_si256(v1), (v0), 1)
#  endif

template <X86_VEXT vext>
static void processBO_SIMD(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset, int startIdx,
                           const ClpRng& clpRng, int channelBitDepth) {
  const int shiftBits = channelBitDepth - NUM_SAO_BO_CLASSES_LOG2;
  int8_t p_eo_offsets[16] = {0};
  for (int i = 0; i < 4; i++) {
    p_eo_offsets[i] = offset[(startIdx + i) % MAX_NUM_SAO_CLASSES];
  }
#  ifdef USE_AVX2
  if (vext >= AVX2) {
    __m256i vbaseoffset = _mm256_set1_epi16(startIdx);
    __m256i vminus = _mm256_set1_epi8(-1);
    __m256i vzero = _mm256_set1_epi8(0);

    __m256i vfour = _mm256_set1_epi16(4);
    __m256i vibdimax = _mm256_set1_epi16((1 << channelBitDepth) - 1);
    __m256i voffsettbl = _mm256_broadcastsi128_si256(_mm_loadu_si128((__m128i*)p_eo_offsets));
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x += 16) {
        __m256i vsrc = _mm256_loadu_si256((__m256i*)&resLine[x]);
        __m256i bands = _mm256_srai_epi16(vsrc, shiftBits);
        bands = _mm256_sub_epi16(bands, vbaseoffset);
        __m256i mask1 = _mm256_cmpgt_epi16(bands, vminus);
        __m256i mask2 = _mm256_cmpgt_epi16(vfour, bands);

        __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, bands);
        veoffsets = _mm256_slli_epi16(veoffsets, 8);
        veoffsets = _mm256_srai_epi16(veoffsets, 8);

        veoffsets = _mm256_and_si256(veoffsets, mask1);
        veoffsets = _mm256_and_si256(veoffsets, mask2);

        vsrc = _mm256_add_epi16(vsrc, veoffsets);
        vsrc = _mm256_min_epi16(_mm256_max_epi16(vsrc, vzero), vibdimax);
        _mm256_storeu_si256((__m256i*)&resLine[x], vsrc);
      }
      resLine += resStride;
    }
    _mm256_zeroupper();
  } else
#  endif
  {
    __m128i vbaseoffset = _mm_set1_epi16(startIdx);
    __m128i vminus = _mm_set1_epi8(-1);
    __m128i vzero = _mm_set1_epi8(0);

    __m128i vfour = _mm_set1_epi16(4);
    __m128i vibdimax = _mm_set1_epi16((1 << channelBitDepth) - 1);
    __m128i voffsettbl = _mm_loadu_si128((__m128i*)p_eo_offsets);
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x += 8) {
        __m128i vsrc = _mm_loadu_si128((__m128i*)&resLine[x]);
        __m128i bands = _mm_srai_epi16(vsrc, shiftBits);
        bands = _mm_sub_epi16(bands, vbaseoffset);
        __m128i mask1 = _mm_cmpgt_epi16(bands, vminus);
        __m128i mask2 = _mm_cmplt_epi16(bands, vfour);

        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, bands);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);

        veoffsets = _mm_and_si128(veoffsets, mask1);
        veoffsets = _mm_and_si128(veoffsets, mask2);

        vsrc = _mm_add_epi16(vsrc, veoffsets);
        vsrc = _mm_min_epi16(_mm_max_epi16(vsrc, vzero), vibdimax);
        _mm_store_si128((__m128i*)&resLine[x], vsrc);
      }
      resLine += resStride;
    }
  }
}

template <X86_VEXT vext>
static void processEO0_SIMD(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                            const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, const Pel* leftLine) {
  int startX = isLeftAvail ? 0 : 1;
  int endX = isRightAvail ? width : (width - 1);
  int x, y;
  int8_t p_eo_offsets[16] = {0};
  for (int i = 0; i < SAO_EO_NUM_CATEGORIES; i++) {
    p_eo_offsets[i] = offset[i];
  }
#  ifdef USE_AVX2
  if (vext >= AVX2) {
    if (isLeftAvail && isRightAvail && !(width & 15)) {
      // fast path
      __m256i vbaseoffset = _mm256_set1_epi16(2);
      __m256i vplusone = _mm256_set1_epi16(1);
      __m256i vzero = _mm256_set1_epi8(0);
      __m256i vibdimax = _mm256_set1_epi16(clpRng.max());
      __m256i voffsettbl = _mm256_broadcastsi128_si256(_mm_loadu_si128((__m128i*)p_eo_offsets));
      for (y = 0; y < height; y++) {
        Pel cm1 = leftLine[y];
        //__m256i pre = _mm256_insert_epi16(vzero, cm1, 0);
        for (x = 0; x < width; x += 16) {
          __m256i vsrca, vsrcal, vsrcar;
          vsrca = _mm256_loadu_si256((__m256i*)&resLine[x]);
          vsrcal = _mm256_loadu_si256((__m256i*)&resLine[x - 1]);
          vsrcar = _mm256_loadu_si256((__m256i*)&resLine[x + 1]);
          // todo: find a better way to insert cm1
          vsrcal = _mm256_insert_epi16(vsrcal, cm1, 0);  // this call is slow!

          vsrcal = _mm256_sub_epi16(vsrca, vsrcal);
          vsrcar = _mm256_sub_epi16(vsrca, vsrcar);
          __m256i vsignl = _mm256_sign_epi16(vplusone, vsrcal);
          __m256i vsignr = _mm256_sign_epi16(vplusone, vsrcar);
          __m256i vsign = _mm256_add_epi16(_mm256_add_epi16(vsignl, vsignr), vbaseoffset);
          __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
          veoffsets = _mm256_slli_epi16(veoffsets, 8);
          veoffsets = _mm256_srai_epi16(veoffsets, 8);
          cm1 = resLine[x + 15];  // load last pixel for next loop
          vsrca = _mm256_add_epi16(vsrca, veoffsets);
          vsrca = _mm256_min_epi16(_mm256_max_epi16(vsrca, vzero), vibdimax);

          _mm256_storeu_si256((__m256i*)&resLine[x], vsrca);
        }
        resLine += resStride;
      }
    } else {
      // normal avx2 path
      __m256i vbaseoffset = _mm256_set1_epi16(2);
      __m256i vplusone = _mm256_set1_epi16(1);
      __m256i vzero = _mm256_set1_epi8(0);
      __m256i vibdimax = _mm256_set1_epi16(clpRng.max());
      __m256i voffsettbl = _mm256_broadcastsi128_si256(_mm_loadu_si128((__m128i*)p_eo_offsets));
      for (y = 0; y < height; y++) {
        Pel cm1;
        if (isLeftAvail)
          cm1 = leftLine[y];
        else
          cm1 = resLine[startX - 1];
        for (x = startX; x < endX - 15; x += 16) {
          __m256i vsrca, vsrcal, vsrcar;
          vsrca = _mm256_loadu_si256((__m256i*)&resLine[x]);
          vsrcal = _mm256_loadu_si256((__m256i*)&resLine[x - 1]);
          vsrcar = _mm256_loadu_si256((__m256i*)&resLine[x + 1]);
          vsrcal = _mm256_insert_epi16(vsrcal, cm1, 0);
          vsrcal = _mm256_sub_epi16(vsrca, vsrcal);
          vsrcar = _mm256_sub_epi16(vsrca, vsrcar);
          __m256i vsignl = _mm256_sign_epi16(vplusone, vsrcal);
          __m256i vsignr = _mm256_sign_epi16(vplusone, vsrcar);
          __m256i vsign = _mm256_add_epi16(_mm256_add_epi16(vsignl, vsignr), vbaseoffset);
          __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
          veoffsets = _mm256_slli_epi16(veoffsets, 8);
          veoffsets = _mm256_srai_epi16(veoffsets, 8);
          cm1 = resLine[x + 15];  // load last pixel for next loop
          vsrca = _mm256_add_epi16(vsrca, veoffsets);
          vsrca = _mm256_min_epi16(_mm256_max_epi16(vsrca, vzero), vibdimax);

          _mm256_storeu_si256((__m256i*)&resLine[x], vsrca);
        }
        if (x < endX) {
          Pel tmp[16];
          __m256i vsrca, vsrcal, vsrcar;
          vsrca = _mm256_loadu_si256((__m256i*)&resLine[x]);
          vsrcal = _mm256_loadu_si256((__m256i*)&resLine[x - 1]);
          vsrcar = _mm256_loadu_si256((__m256i*)&resLine[x + 1]);
          vsrcal = _mm256_insert_epi16(vsrcal, cm1, 0);
          vsrcal = _mm256_sub_epi16(vsrca, vsrcal);
          vsrcar = _mm256_sub_epi16(vsrca, vsrcar);
          __m256i vsignl = _mm256_sign_epi16(vplusone, vsrcal);
          __m256i vsignr = _mm256_sign_epi16(vplusone, vsrcar);
          __m256i vsign = _mm256_add_epi16(_mm256_add_epi16(vsignl, vsignr), vbaseoffset);
          __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
          veoffsets = _mm256_slli_epi16(veoffsets, 8);
          veoffsets = _mm256_srai_epi16(veoffsets, 8);
          vsrca = _mm256_add_epi16(vsrca, veoffsets);
          vsrca = _mm256_min_epi16(_mm256_max_epi16(vsrca, vzero), vibdimax);

          _mm256_storeu_si256((__m256i*)&tmp[0], vsrca);
          int base = x;
          for (; x < endX; x++) {
            resLine[x] = tmp[x - base];
          }
        }
        resLine += resStride;
      }
    }
    _mm256_zeroupper();
  } else
#  endif
  {
    __m128i vsrca, vsrcal, vsrcar;
    __m128i vbaseoffset = _mm_set1_epi16(2);
    __m128i vplusone = _mm_set1_epi16(1);
    __m128i vzero = _mm_set1_epi8(0);
    __m128i vibdimax = _mm_set1_epi16(clpRng.max());
    __m128i voffsettbl = _mm_loadu_si128((__m128i*)p_eo_offsets);
    // sse simd
    for (y = 0; y < height; y++) {
      Pel cm1;
      if (isLeftAvail)
        cm1 = leftLine[y];
      else
        cm1 = resLine[startX - 1];
      for (x = startX; x < endX - 7; x += 8) {
        vsrca = _mm_loadu_si128((__m128i*)&resLine[x]);
        vsrcal = _mm_loadu_si128((__m128i*)&resLine[x - 1]);
        vsrcar = _mm_loadu_si128((__m128i*)&resLine[x + 1]);
        vsrcal = _mm_insert_epi16(vsrcal, cm1, 0);
        vsrcal = _mm_sub_epi16(vsrca, vsrcal);
        vsrcar = _mm_sub_epi16(vsrca, vsrcar);
        __m128i vsignl = _mm_sign_epi16(vplusone, vsrcal);
        __m128i vsignr = _mm_sign_epi16(vplusone, vsrcar);
        __m128i vsign = _mm_add_epi16(_mm_add_epi16(vsignl, vsignr), vbaseoffset);
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        cm1 = resLine[x + 7];
        vsrca = _mm_add_epi16(vsrca, veoffsets);
        vsrca = _mm_min_epi16(_mm_max_epi16(vsrca, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&resLine[x], vsrca);
      }
      if (x < endX) {
        Pel tmp[8];
        int base = x;
        vsrca = _mm_loadu_si128((__m128i*)&resLine[x]);
        vsrcal = _mm_loadu_si128((__m128i*)&resLine[x - 1]);
        vsrcar = _mm_loadu_si128((__m128i*)&resLine[x + 1]);
        vsrcal = _mm_insert_epi16(vsrcal, cm1, 0);
        vsrcal = _mm_sub_epi16(vsrca, vsrcal);
        vsrcar = _mm_sub_epi16(vsrca, vsrcar);
        __m128i vsignl = _mm_sign_epi16(vplusone, vsrcal);
        __m128i vsignr = _mm_sign_epi16(vplusone, vsrcar);
        __m128i vsign = _mm_add_epi16(_mm_add_epi16(vsignl, vsignr), vbaseoffset);
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        vsrca = _mm_add_epi16(vsrca, veoffsets);
        vsrca = _mm_min_epi16(_mm_max_epi16(vsrca, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&tmp[0], vsrca);
        for (; x < endX; x++) {
          resLine[x] = tmp[x - base];
        }
      }
      resLine += resStride;
    }
  }
}

template <X86_VEXT vext>
static void processEO90_SIMD(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                             const ClpRng& clpRng, bool isAboveAvail, bool isBottomAvail, const Pel* topLine) {
  int16_t signDownLine[128];
  int startY = isAboveAvail ? 0 : 1;
  int endY = isBottomAvail ? height : (height - 1);
  int x, y;
  int8_t p_eo_offsets[16] = {0};
  for (int i = 0; i < SAO_EO_NUM_CATEGORIES; i++) {
    p_eo_offsets[i] = offset[i];
  }
  if (!isAboveAvail) {
    resLine += resStride;
  }
  const Pel* srcLineAbove = resLine - resStride;
  if (isAboveAvail) srcLineAbove = topLine;
#  ifdef USE_AVX2
  if (vext >= AVX2) {
    __m256i vbaseoffset = _mm256_set1_epi16(2);
    __m256i vplusone = _mm256_set1_epi16(1);
    __m256i vzero = _mm256_set1_epi8(0);
    __m256i vibdimax = _mm256_set1_epi16(clpRng.max());
    __m256i voffsettbl = _mm256_broadcastsi128_si256(_mm_loadu_si128((__m128i*)p_eo_offsets));
    for (x = 0; x < width; x += 16) {
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i above = _mm256_loadu_si256((__m256i*)&srcLineAbove[x]);
      curr = _mm256_sub_epi16(above, curr);
      curr = _mm256_sign_epi16(vplusone, curr);
      _mm256_storeu_si256((__m256i*)&signDownLine[x], curr);
    }
    for (y = startY; y < endY; y++) {
      const Pel* srcLineBelow = resLine + resStride;
      for (x = 0; x < width; x += 16) {
        __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
        __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x]);
        __m256i signDown = _mm256_loadu_si256((__m256i*)&signDownLine[x]);
        below = _mm256_sub_epi16(curr, below);
        signDown = _mm256_sub_epi16(vbaseoffset, signDown);  // 2 - last signdown
        below = _mm256_sign_epi16(vplusone, below);          // new signDown
        signDown = _mm256_add_epi16(signDown, below);
        __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, signDown);
        veoffsets = _mm256_slli_epi16(veoffsets, 8);
        veoffsets = _mm256_srai_epi16(veoffsets, 8);
        _mm256_storeu_si256((__m256i*)&signDownLine[x], below);
        curr = _mm256_add_epi16(curr, veoffsets);
        curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
        _mm256_storeu_si256((__m256i*)&resLine[x], curr);
      }
      resLine += resStride;
    }
    _mm256_zeroupper();
  } else
#  endif
  {
    __m128i vbaseoffset = _mm_set1_epi16(2);
    __m128i vplusone = _mm_set1_epi16(1);
    __m128i vzero = _mm_set1_epi8(0);
    __m128i vibdimax = _mm_set1_epi16(clpRng.max());
    __m128i voffsettbl = _mm_loadu_si128((__m128i*)p_eo_offsets);
    for (x = 0; x < width; x += 8) {
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i above = _mm_loadu_si128((__m128i*)&srcLineAbove[x]);
      curr = _mm_sub_epi16(above, curr);
      curr = _mm_sign_epi16(vplusone, curr);
      _mm_storeu_si128((__m128i*)&signDownLine[x], curr);
    }
    for (y = startY; y < endY; y++) {
      const Pel* srcLineBelow = resLine + resStride;
      for (x = 0; x < width; x += 8) {
        __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
        __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x]);
        __m128i signDown = _mm_loadu_si128((__m128i*)&signDownLine[x]);
        below = _mm_sub_epi16(curr, below);
        signDown = _mm_sub_epi16(vbaseoffset, signDown);  // 2 - last signdown
        below = _mm_sign_epi16(vplusone, below);          // new signDown
        signDown = _mm_add_epi16(signDown, below);
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, signDown);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        _mm_storeu_si128((__m128i*)&signDownLine[x], below);
        curr = _mm_add_epi16(curr, veoffsets);
        curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&resLine[x], curr);
      }
      resLine += resStride;
    }
  }
}

template <X86_VEXT vext>
static void processEO135_SIMD(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                              const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, bool isAboveLeftAvail,
                              bool isAboveAvail, bool isBelowAvail, bool isBelowRightAvail, const Pel* topLine,
                              const Pel* leftLine) {
  int8_t p_eo_offsets[16] = {0};
  for (int i = 0; i < SAO_EO_NUM_CATEGORIES; i++) {
    p_eo_offsets[i] = offset[i];
  }
  offset += 2;
  int16_t *signUpLine, *signDownLine, *signTmpLine;
  int16_t aSignUpLine[128 + 2], aSignDownLine[128 + 2];
  signUpLine = &aSignUpLine[0];
  signDownLine = &aSignDownLine[0];

  int startX = isLeftAvail ? 0 : 1;
  int endX = isRightAvail ? width : (width - 1);
  int x, y;
  // prepare 2nd line's upper sign
  const Pel* srcLineBelow = resLine + resStride;
  x = startX;
  if (isLeftAvail) {
    signUpLine[x] = (int16_t)sgn(srcLineBelow[x] - leftLine[0]);
    x++;
  }
#  ifdef USE_AVX2
  if (vext >= AVX2) {
    __m256i vbaseoffset = _mm256_set1_epi16(2);
    __m256i vplusone = _mm256_set1_epi16(1);
    __m256i vzero = _mm256_set1_epi8(0);
    __m256i vibdimax = _mm256_set1_epi16(clpRng.max());
    __m256i voffsettbl = _mm256_broadcastsi128_si256(_mm_loadu_si128((__m128i*)p_eo_offsets));

    for (; x < endX + 1; x += 16) {
      __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x]);
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x - 1]);
      below = _mm256_sub_epi16(below, curr);
      below = _mm256_sign_epi16(vplusone, below);
      _mm256_storeu_si256((__m256i*)&signUpLine[x], below);
    }
    // 1st line
    const Pel* srcLineAbove = resLine - resStride;
    int firstLineStartX = isAboveLeftAvail ? 0 : 1;
    int firstLineEndX = isAboveAvail ? endX : 1;
    if (isAboveAvail) srcLineAbove = topLine;
    for (x = firstLineStartX; x < firstLineEndX - 15; x += 16) {
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i above = _mm256_loadu_si256((__m256i*)&srcLineAbove[x - 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x + 1]);
      above = _mm256_sub_epi16(curr, above);
      signUp = _mm256_sub_epi16(vbaseoffset, signUp);
      above = _mm256_sign_epi16(vplusone, above);
      __m256i vsign = _mm256_add_epi16(signUp, above);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);

      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&resLine[x], curr);
    }
    if (x < firstLineEndX) {
      Pel tmpData[16];
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i above = _mm256_loadu_si256((__m256i*)&srcLineAbove[x - 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x + 1]);
      above = _mm256_sub_epi16(curr, above);
      signUp = _mm256_sub_epi16(vbaseoffset, signUp);
      above = _mm256_sign_epi16(vplusone, above);
      __m256i vsign = _mm256_add_epi16(signUp, above);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);

      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&tmpData[0], curr);
      int base = x;
      for (; x < firstLineEndX; x++) resLine[x] = tmpData[x - base];
    }
    resLine += resStride;
    // middle lines
    for (y = 1; y < height - 1; y++) {
      srcLineBelow = resLine + resStride;
      for (x = startX; x < endX - 15; x += 16) {
        __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
        __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x + 1]);
        __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
        below = _mm256_sub_epi16(curr, below);
        signUp = _mm256_add_epi16(signUp, vbaseoffset);
        below = _mm256_sign_epi16(vplusone, below);
        __m256i vsign = _mm256_add_epi16(signUp, below);
        __m256i next = _mm256_sub_epi16(vzero, below);
        __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
        _mm256_storeu_si256((__m256i*)&signDownLine[x + 1], next);
        veoffsets = _mm256_slli_epi16(veoffsets, 8);
        veoffsets = _mm256_srai_epi16(veoffsets, 8);
        curr = _mm256_add_epi16(curr, veoffsets);
        curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
        _mm256_storeu_si256((__m256i*)&resLine[x], curr);
      }
      if (x < endX) {
        Pel tmpData[16];
        __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
        __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x + 1]);
        __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
        below = _mm256_sub_epi16(curr, below);
        signUp = _mm256_add_epi16(signUp, vbaseoffset);
        below = _mm256_sign_epi16(vplusone, below);
        __m256i vsign = _mm256_add_epi16(signUp, below);
        __m256i next = _mm256_sub_epi16(vzero, below);
        __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
        _mm256_storeu_si256((__m256i*)&signDownLine[x + 1], next);
        veoffsets = _mm256_slli_epi16(veoffsets, 8);
        veoffsets = _mm256_srai_epi16(veoffsets, 8);
        curr = _mm256_add_epi16(curr, veoffsets);
        curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
        _mm256_storeu_si256((__m256i*)&tmpData[0], curr);
        int base = x;
        for (; x < endX; x++) resLine[x] = tmpData[x - base];
      }
      if (isLeftAvail)
        signDownLine[startX] = (int16_t)sgn(srcLineBelow[startX] - leftLine[y]);
      else
        signDownLine[startX] = (int16_t)sgn(srcLineBelow[startX] - resLine[startX - 1]);
      signTmpLine = signUpLine;
      signUpLine = signDownLine;
      signDownLine = signTmpLine;
      resLine += resStride;
    }

    // last line
    srcLineBelow = resLine + resStride;
    int lastLineStartX = isBelowAvail ? startX : (width - 1);
    int lastLineEndX = isBelowRightAvail ? width : (width - 1);
    for (x = lastLineStartX; x < lastLineEndX - 15; x += 16) {
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x + 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
      below = _mm256_sub_epi16(curr, below);
      signUp = _mm256_add_epi16(signUp, vbaseoffset);
      below = _mm256_sign_epi16(vplusone, below);
      __m256i vsign = _mm256_add_epi16(signUp, below);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);
      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&resLine[x], curr);
    }
    if (x < lastLineEndX) {
      Pel tmpData[16];
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x + 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
      below = _mm256_sub_epi16(curr, below);
      signUp = _mm256_add_epi16(signUp, vbaseoffset);
      below = _mm256_sign_epi16(vplusone, below);
      __m256i vsign = _mm256_add_epi16(signUp, below);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);
      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&tmpData[0], curr);
      int base = x;
      for (; x < lastLineEndX; x++) resLine[x] = tmpData[x - base];
    }
  } else
#  endif
  {
    __m128i vbaseoffset = _mm_set1_epi16(2);
    __m128i vplusone = _mm_set1_epi16(1);
    __m128i vzero = _mm_set1_epi8(0);
    __m128i vibdimax = _mm_set1_epi16(clpRng.max());
    __m128i voffsettbl = _mm_loadu_si128((__m128i*)p_eo_offsets);

    for (; x < endX + 1; x += 8) {
      __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x]);
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x - 1]);
      below = _mm_sub_epi16(below, curr);
      below = _mm_sign_epi16(vplusone, below);
      _mm_storeu_si128((__m128i*)&signUpLine[x], below);
    }
    // 1st line
    const Pel* srcLineAbove = resLine - resStride;
    int firstLineStartX = isAboveLeftAvail ? 0 : 1;
    int firstLineEndX = isAboveAvail ? endX : 1;
    if (isAboveAvail) srcLineAbove = topLine;
    for (x = firstLineStartX; x < firstLineEndX - 7; x += 8) {
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i above = _mm_loadu_si128((__m128i*)&srcLineAbove[x - 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x + 1]);
      above = _mm_sub_epi16(curr, above);
      signUp = _mm_sub_epi16(vbaseoffset, signUp);
      above = _mm_sign_epi16(vplusone, above);
      __m128i vsign = _mm_add_epi16(signUp, above);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);

      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&resLine[x], curr);
    }
    if (x < firstLineEndX) {
      Pel tmpData[8];
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i above = _mm_loadu_si128((__m128i*)&srcLineAbove[x - 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x + 1]);
      above = _mm_sub_epi16(curr, above);
      signUp = _mm_sub_epi16(vbaseoffset, signUp);
      above = _mm_sign_epi16(vplusone, above);
      __m128i vsign = _mm_add_epi16(signUp, above);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);

      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&tmpData[0], curr);
      int base = x;
      for (; x < firstLineEndX; x++) resLine[x] = tmpData[x - base];
    }
    resLine += resStride;
    // middle lines
    for (y = 1; y < height - 1; y++) {
      srcLineBelow = resLine + resStride;
      for (x = startX; x < endX - 7; x += 8) {
        __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
        __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x + 1]);
        __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
        below = _mm_sub_epi16(curr, below);
        signUp = _mm_add_epi16(signUp, vbaseoffset);
        below = _mm_sign_epi16(vplusone, below);
        __m128i vsign = _mm_add_epi16(signUp, below);
        __m128i next = _mm_sub_epi16(vzero, below);
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
        _mm_storeu_si128((__m128i*)&signDownLine[x + 1], next);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        curr = _mm_add_epi16(curr, veoffsets);
        curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&resLine[x], curr);
      }
      if (x < endX) {
        Pel tmpData[8];
        __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
        __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x + 1]);
        __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
        below = _mm_sub_epi16(curr, below);
        signUp = _mm_add_epi16(signUp, vbaseoffset);
        below = _mm_sign_epi16(vplusone, below);
        __m128i vsign = _mm_add_epi16(signUp, below);
        __m128i next = _mm_sub_epi16(vzero, below);
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
        _mm_storeu_si128((__m128i*)&signDownLine[x + 1], next);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        curr = _mm_add_epi16(curr, veoffsets);
        curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&tmpData[0], curr);
        int base = x;
        for (; x < endX; x++) resLine[x] = tmpData[x - base];
      }
      if (isLeftAvail)
        signDownLine[startX] = (int16_t)sgn(srcLineBelow[startX] - leftLine[y]);
      else
        signDownLine[startX] = (int16_t)sgn(srcLineBelow[startX] - resLine[startX - 1]);
      signTmpLine = signUpLine;
      signUpLine = signDownLine;
      signDownLine = signTmpLine;
      resLine += resStride;
    }

    // last line
    srcLineBelow = resLine + resStride;
    int lastLineStartX = isBelowAvail ? startX : (width - 1);
    int lastLineEndX = isBelowRightAvail ? width : (width - 1);
    for (x = lastLineStartX; x < lastLineEndX - 7; x += 8) {
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x + 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
      below = _mm_sub_epi16(curr, below);
      signUp = _mm_add_epi16(signUp, vbaseoffset);
      below = _mm_sign_epi16(vplusone, below);
      __m128i vsign = _mm_add_epi16(signUp, below);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);
      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&resLine[x], curr);
    }
    if (x < lastLineEndX) {
      Pel tmpData[8];
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x + 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
      below = _mm_sub_epi16(curr, below);
      signUp = _mm_add_epi16(signUp, vbaseoffset);
      below = _mm_sign_epi16(vplusone, below);
      __m128i vsign = _mm_add_epi16(signUp, below);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);
      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&tmpData[0], curr);
      int base = x;
      for (; x < lastLineEndX; x++) resLine[x] = tmpData[x - base];
    }
  }
}

template <X86_VEXT vext>
static void processEO45_SIMD(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                             const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, bool isAboveRightAvail,
                             bool isAboveAvail, bool isBelowAvail, bool isBelowLeftAvail, const Pel* topLine,
                             const Pel* leftLine) {
  int8_t p_eo_offsets[16] = {0};
  for (int i = 0; i < SAO_EO_NUM_CATEGORIES; i++) {
    p_eo_offsets[i] = offset[i];
  }
  offset += 2;
  int16_t aSignUpLine[128 + 2];
  int16_t* signUpLine = &aSignUpLine[1];
  int x, y;
  int startX = isLeftAvail ? 0 : 1;
  int endX = isRightAvail ? width : (width - 1);

  // prepare 2nd line upper sign
  const Pel* srcLineBelow = resLine + resStride;
  x = startX - 1;
  if (isLeftAvail) {
    signUpLine[x] = (int16_t)sgn(leftLine[1] - resLine[x + 1]);
    x++;
  }
#  ifdef USE_AVX2
  if (vext >= AVX2) {
    __m256i vbaseoffset = _mm256_set1_epi16(2);
    __m256i vplusone = _mm256_set1_epi16(1);
    __m256i vzero = _mm256_set1_epi8(0);
    __m256i vibdimax = _mm256_set1_epi16(clpRng.max());
    __m256i voffsettbl = _mm256_broadcastsi128_si256(_mm_loadu_si128((__m128i*)p_eo_offsets));
    for (; x < endX; x += 16) {
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x + 1]);
      __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x]);
      below = _mm256_sub_epi16(below, curr);
      below = _mm256_sign_epi16(vplusone, below);
      _mm256_storeu_si256((__m256i*)&signUpLine[x], below);
    }
    // for (; x< endX; x++)
    //{
    //  signUpLine[x] = (int16_t)sgn(srcLineBelow[x] - resLine[x+1]);
    //}
    // first line
    const Pel* srcLineAbove = resLine - resStride;
    int firstLineStartX = isAboveAvail ? startX : (width - 1);
    int firstLineEndX = isAboveRightAvail ? width : (width - 1);
    if (isAboveAvail) srcLineAbove = topLine;
    for (x = firstLineStartX; x < firstLineEndX - 15; x += 16) {
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i above = _mm256_loadu_si256((__m256i*)&srcLineAbove[x + 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x - 1]);
      above = _mm256_sub_epi16(curr, above);
      signUp = _mm256_sub_epi16(vbaseoffset, signUp);  // 2 - signUpLine[x-1];
      above = _mm256_sign_epi16(vplusone, above);      // sgn(resLine[x] - srcLineAbove[x+1])
      __m256i vsign = _mm256_add_epi16(above, signUp);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);
      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&resLine[x], curr);
    }
    if (x < firstLineEndX) {
      Pel tmpData[16];
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i above = _mm256_loadu_si256((__m256i*)&srcLineAbove[x + 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x - 1]);
      above = _mm256_sub_epi16(curr, above);
      signUp = _mm256_sub_epi16(vbaseoffset, signUp);  // 2 - signUpLine[x-1];
      above = _mm256_sign_epi16(vplusone, above);      // sgn(resLine[x] - srcLineAbove[x+1])
      __m256i vsign = _mm256_add_epi16(above, signUp);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);
      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&tmpData[0], curr);
      int base = x;
      for (; x < firstLineEndX; x++) resLine[x] = tmpData[x - base];
    }
    // for(x= firstLineStartX; x< firstLineEndX; x++)
    //{
    //  int edgeType = sgn(resLine[x] - srcLineAbove[x+1]) - signUpLine[x-1];
    //  resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
    //}
    resLine += resStride;

    // middle lines
    for (y = 1; y < height - 1; y++) {
      srcLineBelow = resLine + resStride;
      x = startX;
      if (isLeftAvail) {
        int signDown = (int16_t)sgn(resLine[x] - leftLine[y + 1]);
        int edgeType = signDown + signUpLine[x];
        resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
        signUpLine[x - 1] = -signDown;
        x++;
      }
      for (; x < endX - 15; x += 16) {
        __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
        __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x - 1]);
        __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
        __m256i signDown = _mm256_sub_epi16(curr, below);
        signUp = _mm256_add_epi16(signUp, vbaseoffset);
        signDown = _mm256_sign_epi16(vplusone, signDown);
        __m256i vsign = _mm256_add_epi16(signUp, signDown);
        below = _mm256_sub_epi16(vzero, signDown);  // -signDown
        __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
        _mm256_storeu_si256((__m256i*)&signUpLine[x - 1], below);
        veoffsets = _mm256_slli_epi16(veoffsets, 8);
        veoffsets = _mm256_srai_epi16(veoffsets, 8);
        curr = _mm256_add_epi16(curr, veoffsets);
        curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
        _mm256_storeu_si256((__m256i*)&resLine[x], curr);
      }
      if (x < endX) {
        Pel tmpData[16];
        __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
        __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x - 1]);
        __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
        __m256i signDown = _mm256_sub_epi16(curr, below);
        signUp = _mm256_add_epi16(signUp, vbaseoffset);
        signDown = _mm256_sign_epi16(vplusone, signDown);
        __m256i vsign = _mm256_add_epi16(signUp, signDown);
        below = _mm256_sub_epi16(vzero, signDown);  // -signDown
        __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
        _mm256_storeu_si256((__m256i*)&signUpLine[x - 1], below);
        veoffsets = _mm256_slli_epi16(veoffsets, 8);
        veoffsets = _mm256_srai_epi16(veoffsets, 8);
        curr = _mm256_add_epi16(curr, veoffsets);
        curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
        _mm256_storeu_si256((__m256i*)&tmpData[0], curr);
        int base = x;
        for (; x < endX; x++) resLine[x] = tmpData[x - base];
      }
      // for(; x< endX; x++)
      //{
      //  int signDown =  (int16_t)sgn(resLine[x] - srcLineBelow[x-1]);
      //  int edgeType =  signDown + signUpLine[x];
      //  resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      //  signUpLine[x-1] = -signDown;
      //}
      signUpLine[endX - 1] = (int16_t)sgn(srcLineBelow[endX - 1] - resLine[endX]);
      resLine += resStride;
    }
    // last line
    srcLineBelow = resLine + resStride;
    int lastLineStartX = isBelowLeftAvail ? 0 : 1;
    int lastLineEndX = isBelowAvail ? endX : 1;
    // below left is in reconstruction buffer, not line buffer
    for (x = lastLineStartX; x < lastLineEndX - 15; x += 16) {
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x - 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
      __m256i signDown = _mm256_sub_epi16(curr, below);
      signUp = _mm256_add_epi16(signUp, vbaseoffset);
      signDown = _mm256_sign_epi16(vplusone, signDown);
      __m256i vsign = _mm256_add_epi16(signUp, signDown);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);
      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&resLine[x], curr);
    }
    if (x < lastLineEndX) {
      Pel tmpData[16];
      __m256i curr = _mm256_loadu_si256((__m256i*)&resLine[x]);
      __m256i below = _mm256_loadu_si256((__m256i*)&srcLineBelow[x - 1]);
      __m256i signUp = _mm256_loadu_si256((__m256i*)&signUpLine[x]);
      __m256i signDown = _mm256_sub_epi16(curr, below);
      signUp = _mm256_add_epi16(signUp, vbaseoffset);
      signDown = _mm256_sign_epi16(vplusone, signDown);
      __m256i vsign = _mm256_add_epi16(signUp, signDown);
      __m256i veoffsets = _mm256_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm256_slli_epi16(veoffsets, 8);
      veoffsets = _mm256_srai_epi16(veoffsets, 8);
      curr = _mm256_add_epi16(curr, veoffsets);
      curr = _mm256_min_epi16(_mm256_max_epi16(curr, vzero), vibdimax);
      _mm256_storeu_si256((__m256i*)&tmpData[0], curr);
      int base = x;
      for (; x < lastLineEndX; x++) resLine[x] = tmpData[x - base];
    }
    // for(x = lastLineStartX; x< lastLineEndX; x++)
    //{
    //  int edgeType = sgn(resLine[x] - srcLineBelow[x-1]) + signUpLine[x];
    //  resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
    //}
  } else
#  endif
  {
    __m128i vbaseoffset = _mm_set1_epi16(2);
    __m128i vplusone = _mm_set1_epi16(1);
    __m128i vzero = _mm_set1_epi8(0);
    __m128i vibdimax = _mm_set1_epi16(clpRng.max());
    __m128i voffsettbl = _mm_loadu_si128((__m128i*)p_eo_offsets);
    for (; x < endX; x += 8) {
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x + 1]);
      __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x]);
      below = _mm_sub_epi16(below, curr);
      below = _mm_sign_epi16(vplusone, below);
      _mm_storeu_si128((__m128i*)&signUpLine[x], below);
    }
    // first line
    const Pel* srcLineAbove = resLine - resStride;
    int firstLineStartX = isAboveAvail ? startX : (width - 1);
    int firstLineEndX = isAboveRightAvail ? width : (width - 1);
    if (isAboveAvail) srcLineAbove = topLine;
    for (x = firstLineStartX; x < firstLineEndX - 7; x += 8) {
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i above = _mm_loadu_si128((__m128i*)&srcLineAbove[x + 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x - 1]);
      above = _mm_sub_epi16(curr, above);
      signUp = _mm_sub_epi16(vbaseoffset, signUp);  // 2 - signUpLine[x-1];
      above = _mm_sign_epi16(vplusone, above);      // sgn(resLine[x] - srcLineAbove[x+1])
      __m128i vsign = _mm_add_epi16(above, signUp);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);
      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&resLine[x], curr);
    }
    if (x < firstLineEndX) {
      Pel tmpData[8];
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i above = _mm_loadu_si128((__m128i*)&srcLineAbove[x + 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x - 1]);
      above = _mm_sub_epi16(curr, above);
      signUp = _mm_sub_epi16(vbaseoffset, signUp);  // 2 - signUpLine[x-1];
      above = _mm_sign_epi16(vplusone, above);      // sgn(resLine[x] - srcLineAbove[x+1])
      __m128i vsign = _mm_add_epi16(above, signUp);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);
      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&tmpData[0], curr);
      int base = x;
      for (; x < firstLineEndX; x++) resLine[x] = tmpData[x - base];
    }
    resLine += resStride;

    // middle lines
    for (y = 1; y < height - 1; y++) {
      srcLineBelow = resLine + resStride;
      x = startX;
      if (isLeftAvail) {
        int signDown = (int16_t)sgn(resLine[x] - leftLine[y + 1]);
        int edgeType = signDown + signUpLine[x];
        resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
        signUpLine[x - 1] = -signDown;
        x++;
      }
      for (; x < endX - 7; x += 8) {
        __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
        __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x - 1]);
        __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
        __m128i signDown = _mm_sub_epi16(curr, below);
        signUp = _mm_add_epi16(signUp, vbaseoffset);
        signDown = _mm_sign_epi16(vplusone, signDown);
        __m128i vsign = _mm_add_epi16(signUp, signDown);
        below = _mm_sub_epi16(vzero, signDown);  // -signDown
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
        _mm_storeu_si128((__m128i*)&signUpLine[x - 1], below);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        curr = _mm_add_epi16(curr, veoffsets);
        curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&resLine[x], curr);
      }
      if (x < endX) {
        Pel tmpData[8];
        __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
        __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x - 1]);
        __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
        __m128i signDown = _mm_sub_epi16(curr, below);
        signUp = _mm_add_epi16(signUp, vbaseoffset);
        signDown = _mm_sign_epi16(vplusone, signDown);
        __m128i vsign = _mm_add_epi16(signUp, signDown);
        below = _mm_sub_epi16(vzero, signDown);  // -signDown
        __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
        _mm_storeu_si128((__m128i*)&signUpLine[x - 1], below);
        veoffsets = _mm_slli_epi16(veoffsets, 8);
        veoffsets = _mm_srai_epi16(veoffsets, 8);
        curr = _mm_add_epi16(curr, veoffsets);
        curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
        _mm_storeu_si128((__m128i*)&tmpData[0], curr);
        int base = x;
        for (; x < endX; x++) resLine[x] = tmpData[x - base];
      }
      signUpLine[endX - 1] = (int16_t)sgn(srcLineBelow[endX - 1] - resLine[endX]);
      resLine += resStride;
    }
    // last line
    srcLineBelow = resLine + resStride;
    int lastLineStartX = isBelowLeftAvail ? 0 : 1;
    int lastLineEndX = isBelowAvail ? endX : 1;
    // below left is in reconstruction buffer, not line buffer
    for (x = lastLineStartX; x < lastLineEndX - 7; x += 8) {
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x - 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
      __m128i signDown = _mm_sub_epi16(curr, below);
      signUp = _mm_add_epi16(signUp, vbaseoffset);
      signDown = _mm_sign_epi16(vplusone, signDown);
      __m128i vsign = _mm_add_epi16(signUp, signDown);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);
      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&resLine[x], curr);
    }
    if (x < lastLineEndX) {
      Pel tmpData[8];
      __m128i curr = _mm_loadu_si128((__m128i*)&resLine[x]);
      __m128i below = _mm_loadu_si128((__m128i*)&srcLineBelow[x - 1]);
      __m128i signUp = _mm_loadu_si128((__m128i*)&signUpLine[x]);
      __m128i signDown = _mm_sub_epi16(curr, below);
      signUp = _mm_add_epi16(signUp, vbaseoffset);
      signDown = _mm_sign_epi16(vplusone, signDown);
      __m128i vsign = _mm_add_epi16(signUp, signDown);
      __m128i veoffsets = _mm_shuffle_epi8(voffsettbl, vsign);
      veoffsets = _mm_slli_epi16(veoffsets, 8);
      veoffsets = _mm_srai_epi16(veoffsets, 8);
      curr = _mm_add_epi16(curr, veoffsets);
      curr = _mm_min_epi16(_mm_max_epi16(curr, vzero), vibdimax);
      _mm_storeu_si128((__m128i*)&tmpData[0], curr);
      int base = x;
      for (; x < lastLineEndX; x++) resLine[x] = tmpData[x - base];
    }
  }
}

template <X86_VEXT vext>
void SampleAdaptiveOffset::_initSampleAdaptiveOffsetX86(int bitDepth) {
  processBO = processBO_SIMD<vext>;
  processEO0 = processEO0_SIMD<vext>;
  processEO90 = processEO90_SIMD<vext>;
  processEO135 = processEO135_SIMD<vext>;
  processEO45 = processEO45_SIMD<vext>;
}

template void SampleAdaptiveOffset::_initSampleAdaptiveOffsetX86<SIMDX86>(int bitDepth);
#endif  //#ifdef TARGET_SIMD_X86
//! \}
