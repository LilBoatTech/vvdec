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

/** \file     RdCostX86.cpp
    \brief    RD cost computation class, SIMD version
*/

#include <math.h>
#include <limits>

#include "CommonDefX86.h"
#include "../RdCost.h"

#ifdef TARGET_SIMD_X86

template <X86_VEXT vext>
Distortion RdCost::xGetSAD_16xN_SIMD(const DistParam& rcDtParam) {
  if (rcDtParam.bitDepth > 10) return RdCost::xGetSAD(rcDtParam);

  //  assert( rcDtParam.iCols == iWidth);
  const short* pSrc1 = (const short*)rcDtParam.org.buf;
  const short* pSrc2 = (const short*)rcDtParam.cur.buf;
  const int iRows = rcDtParam.org.height;
  const int iSubShift = rcDtParam.subShift;
  const ptrdiff_t iStrideSrc1 = rcDtParam.org.stride << iSubShift;
  const ptrdiff_t iStrideSrc2 = rcDtParam.cur.stride << iSubShift;

  uint32_t uiSum = 0;

  if (vext >= AVX2) {
#  ifdef USE_AVX2
    __m256i vone = _mm256_set1_epi16(1);
    __m256i vsum32 = _mm256_setzero_si256();

    if (iRows == 16 && iSubShift == 1) {
      for (int i = 0; i < 2; i++) {
        __m256i vsum16 = _mm256_setzero_si256();

        // 0
        __m256i vsrc1 = _mm256_loadu_si256((__m256i*)(pSrc1));
        __m256i vsrc2 = _mm256_loadu_si256((__m256i*)(pSrc2));

        pSrc1 += iStrideSrc1;
        pSrc2 += iStrideSrc2;

        // 2x 12bit abs-diff -> 12 bit
        vsum16 = _mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2));

        // 1
        vsrc1 = _mm256_loadu_si256((__m256i*)(pSrc1));
        vsrc2 = _mm256_loadu_si256((__m256i*)(pSrc2));
        pSrc1 += iStrideSrc1;
        pSrc2 += iStrideSrc2;
        // 2x 12bit sum -> 13 bit
        vsum16 = _mm256_add_epi16(vsum16, _mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2)));

        // 2
        vsrc1 = _mm256_loadu_si256((__m256i*)(pSrc1));
        vsrc2 = _mm256_loadu_si256((__m256i*)(pSrc2));
        pSrc1 += iStrideSrc1;
        pSrc2 += iStrideSrc2;
        // 13 bit and 12bit sum, or 3x 12 bit sum -> 13 bit
        vsum16 = _mm256_add_epi16(vsum16, _mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2)));

        // 3
        vsrc1 = _mm256_loadu_si256((__m256i*)(pSrc1));
        vsrc2 = _mm256_loadu_si256((__m256i*)(pSrc2));
        pSrc1 += iStrideSrc1;
        pSrc2 += iStrideSrc2;
        // 4x 12 bit sum -> 14 bit
        vsum16 = _mm256_add_epi16(vsum16, _mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2)));

        vsum32 = _mm256_add_epi32(vsum32, _mm256_madd_epi16(vsum16, vone));
      }
    } else {
      // Do for width that multiple of 16
      for (int iY = 0; iY < iRows; iY += (1 << iSubShift)) {
        __m256i vsrc1 = _mm256_loadu_si256((__m256i*)(pSrc1));
        __m256i vsrc2 = _mm256_loadu_si256((__m256i*)(pSrc2));

        vsum32 = _mm256_add_epi32(vsum32, _mm256_madd_epi16(_mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2)), vone));

        pSrc1 += iStrideSrc1;
        pSrc2 += iStrideSrc2;
      }
    }
    vsum32 = _mm256_hadd_epi32(vsum32, vone);
    vsum32 = _mm256_hadd_epi32(vsum32, vone);
    uiSum = _mm_cvtsi128_si32(_mm256_castsi256_si128(vsum32)) + _mm_cvtsi128_si32(_mm256_extracti128_si256(vsum32, 1));
#  endif
  } else {
    // For width that multiple of 8
    __m128i vone = _mm_set1_epi16(1);
    __m128i vsum32 = _mm_setzero_si128();
    for (int iY = 0; iY < iRows; iY += (1 << iSubShift)) {
      __m128i vsum16 = _mm_setzero_si128();
      for (int iX = 0; iX < 16; iX += 8) {
        __m128i vsrc1 = _mm_loadu_si128((const __m128i*)(&pSrc1[iX]));
        __m128i vsrc2 = _mm_lddqu_si128((const __m128i*)(&pSrc2[iX]));
        vsum16 = _mm_add_epi16(vsum16, _mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2)));
      }
      __m128i vsumtemp = _mm_madd_epi16(vsum16, vone);
      vsum32 = _mm_add_epi32(vsum32, vsumtemp);
      pSrc1 += iStrideSrc1;
      pSrc2 += iStrideSrc2;
    }
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    uiSum = _mm_cvtsi128_si32(vsum32);
  }

  uiSum <<= iSubShift;
  return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth);
}

template <X86_VEXT vext>
Distortion RdCost::xGetSAD_4xN_SIMD(const DistParam& rcDtParam) {
  //  assert( rcDtParam.iCols == iWidth);
  const short* pSrc1 = (const short*)rcDtParam.org.buf;
  const short* pSrc2 = (const short*)rcDtParam.cur.buf;
  int iRows = rcDtParam.org.height;
  int iSubShift = rcDtParam.subShift;
  int iSubStep = (1 << iSubShift);
  const ptrdiff_t iStrideSrc1 = rcDtParam.org.stride * iSubStep;
  const ptrdiff_t iStrideSrc2 = rcDtParam.cur.stride * iSubStep;

  uint32_t uiSum = 0;

  if (iRows == 4 && iSubShift == 0) {
    __m128i vzero = _mm_setzero_si128();
    __m128i vsrc1 = _mm_or_si128(_mm_loadl_epi64((const __m128i*)pSrc1),
                                 _mm_slli_si128(_mm_loadl_epi64((const __m128i*)(&pSrc1[iStrideSrc1])), 8));
    __m128i vsrc2 = _mm_or_si128(_mm_loadl_epi64((const __m128i*)pSrc2),
                                 _mm_slli_si128(_mm_loadl_epi64((const __m128i*)(&pSrc2[iStrideSrc2])), 8));
    __m128i vsum = _mm_cvtepi16_epi32(_mm_hadd_epi16(_mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2)), vzero));

    vsrc1 = _mm_or_si128(_mm_loadl_epi64((const __m128i*)(&pSrc1[2 * iStrideSrc1])),
                         _mm_slli_si128(_mm_loadl_epi64((const __m128i*)(&pSrc1[3 * iStrideSrc1])), 8));
    vsrc2 = _mm_or_si128(_mm_loadl_epi64((const __m128i*)(&pSrc2[2 * iStrideSrc2])),
                         _mm_slli_si128(_mm_loadl_epi64((const __m128i*)(&pSrc2[3 * iStrideSrc2])), 8));
    vsum = _mm_add_epi32(vsum, _mm_cvtepi16_epi32(_mm_hadd_epi16(_mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2)), vzero)));
    vsum = _mm_hadd_epi32(vsum, vzero);
    vsum = _mm_hadd_epi32(vsum, vzero);

    uiSum = _mm_cvtsi128_si32(vsum);
  } else {
    __m128i vone = _mm_set1_epi16(1);
    __m128i vsum32 = _mm_setzero_si128();
    for (int iY = 0; iY < iRows; iY += iSubStep) {
      __m128i vsrc1 = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)pSrc1));
      __m128i vsrc2 = _mm_cvtepi16_epi32(_mm_loadl_epi64((const __m128i*)pSrc2));
      vsum32 = _mm_add_epi32(vsum32, _mm_abs_epi32(_mm_sub_epi32(vsrc1, vsrc2)));

      pSrc1 += iStrideSrc1;
      pSrc2 += iStrideSrc2;
    }
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    uiSum = _mm_cvtsi128_si32(vsum32);
  }

  uiSum <<= iSubShift;
  return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth);
}

template <X86_VEXT vext>
Distortion RdCost::xGetSAD_8xN_SIMD(const DistParam& rcDtParam) {
  //  assert( rcDtParam.iCols == iWidth);
  const short* pSrc1 = (const short*)rcDtParam.org.buf;
  const short* pSrc2 = (const short*)rcDtParam.cur.buf;
  int iRows = rcDtParam.org.height;
  int iSubShift = rcDtParam.subShift;

  if (iRows == 8 && iSubShift == 0) {
    const ptrdiff_t iStrideSrc1 = rcDtParam.org.stride;
    const ptrdiff_t iStrideSrc2 = rcDtParam.cur.stride;

    __m128i vsrc0 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc1 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    __m128i vsrc2 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc3 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    __m128i vsrc4 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc5 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    __m128i vsrc6 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc7 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;

    __m128i vdiff0 = _mm_sub_epi16(vsrc0, vsrc1);
    __m128i vdiff1 = _mm_sub_epi16(vsrc2, vsrc3);
    __m128i vdiff2 = _mm_sub_epi16(vsrc4, vsrc5);
    __m128i vdiff3 = _mm_sub_epi16(vsrc6, vsrc7);

    vdiff0 = _mm_abs_epi16(vdiff0);
    vdiff1 = _mm_abs_epi16(vdiff1);
    vdiff2 = _mm_abs_epi16(vdiff2);
    vdiff3 = _mm_abs_epi16(vdiff3);

    vdiff0 = _mm_add_epi16(vdiff0, vdiff1);
    vdiff2 = _mm_add_epi16(vdiff2, vdiff3);

    __m128i vsum0 = _mm_add_epi16(vdiff0, vdiff2);

    vsrc0 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    vsrc1 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    vsrc2 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    vsrc3 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    vsrc4 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    vsrc5 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    vsrc6 = _mm_loadu_si128((const __m128i*)(pSrc1));
    vsrc7 = _mm_loadu_si128((const __m128i*)(pSrc2));

    vdiff0 = _mm_sub_epi16(vsrc0, vsrc1);
    vdiff1 = _mm_sub_epi16(vsrc2, vsrc3);
    vdiff2 = _mm_sub_epi16(vsrc4, vsrc5);
    vdiff3 = _mm_sub_epi16(vsrc6, vsrc7);

    vdiff0 = _mm_abs_epi16(vdiff0);
    vdiff1 = _mm_abs_epi16(vdiff1);
    vdiff2 = _mm_abs_epi16(vdiff2);
    vdiff3 = _mm_abs_epi16(vdiff3);

    vdiff0 = _mm_add_epi16(vdiff0, vdiff1);
    vdiff2 = _mm_add_epi16(vdiff2, vdiff3);

    __m128i vsum1 = _mm_add_epi16(vdiff0, vdiff2);

    __m128i vone = _mm_set1_epi16(1);
    vsum0 = _mm_add_epi16(vsum0, vsum1);
    __m128i vsum32 = _mm_madd_epi16(vsum0, vone);

    vsum32 = _mm_hadd_epi32(vsum32, vone);
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    uint32_t uiSum = _mm_cvtsi128_si32(vsum32);

    uiSum <<= iSubShift;
    return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth);
  } else if (iRows == 8 && iSubShift == 1) {
    const ptrdiff_t iStrideSrc1 = rcDtParam.org.stride << 1;
    const ptrdiff_t iStrideSrc2 = rcDtParam.cur.stride << 1;

    __m128i vsrc0 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc1 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    __m128i vsrc2 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc3 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    __m128i vsrc4 = _mm_loadu_si128((const __m128i*)(pSrc1));
    pSrc1 += iStrideSrc1;
    __m128i vsrc5 = _mm_loadu_si128((const __m128i*)(pSrc2));
    pSrc2 += iStrideSrc2;
    __m128i vsrc6 = _mm_loadu_si128((const __m128i*)(pSrc1));
    __m128i vsrc7 = _mm_loadu_si128((const __m128i*)(pSrc2));

    __m128i vdiff0 = _mm_sub_epi16(vsrc0, vsrc1);
    __m128i vdiff1 = _mm_sub_epi16(vsrc2, vsrc3);
    __m128i vdiff2 = _mm_sub_epi16(vsrc4, vsrc5);
    __m128i vdiff3 = _mm_sub_epi16(vsrc6, vsrc7);

    vdiff0 = _mm_abs_epi16(vdiff0);
    vdiff1 = _mm_abs_epi16(vdiff1);
    vdiff2 = _mm_abs_epi16(vdiff2);
    vdiff3 = _mm_abs_epi16(vdiff3);

    vdiff0 = _mm_add_epi16(vdiff0, vdiff1);
    vdiff2 = _mm_add_epi16(vdiff2, vdiff3);

    __m128i vone = _mm_set1_epi16(1);
    vdiff0 = _mm_add_epi16(vdiff0, vdiff2);
    __m128i vsum32 = _mm_madd_epi16(vdiff0, vone);

    vsum32 = _mm_hadd_epi32(vsum32, vone);
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    uint32_t uiSum = _mm_cvtsi128_si32(vsum32);

    uiSum <<= iSubShift;
    return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth);
  } else {
    int iSubStep = (1 << iSubShift);
    const ptrdiff_t iStrideSrc1 = rcDtParam.org.stride * iSubStep;
    const ptrdiff_t iStrideSrc2 = rcDtParam.cur.stride * iSubStep;

    uint32_t uiSum = 0;

    // For width that multiple of 8
    __m128i vone = _mm_set1_epi16(1);
    __m128i vsum32 = _mm_setzero_si128();

    for (int iY = 0; iY < iRows; iY += iSubStep) {
      __m128i vsrc1 = _mm_loadu_si128((const __m128i*)(pSrc1));
      __m128i vsrc2 = _mm_loadu_si128((const __m128i*)(pSrc2));
      __m128i vsum16 = _mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2));

      __m128i vsumtemp = _mm_madd_epi16(vsum16, vone);
      vsum32 = _mm_add_epi32(vsum32, vsumtemp);

      pSrc1 += iStrideSrc1;
      pSrc2 += iStrideSrc2;
    }
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    uiSum = _mm_cvtsi128_si32(vsum32);

    uiSum <<= iSubShift;
    return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth);
  }
}

template <int iWidth, X86_VEXT vext>
Distortion RdCost::xGetSAD_16NxN_SIMD(const DistParam& rcDtParam) {
  //  assert( rcDtParam.iCols == iWidth);
  const short* pSrc1 = (const short*)rcDtParam.org.buf;
  const short* pSrc2 = (const short*)rcDtParam.cur.buf;
  int iRows = rcDtParam.org.height;
  int iSubShift = rcDtParam.subShift;
  int iSubStep = (1 << iSubShift);
  const ptrdiff_t iStrideSrc1 = rcDtParam.org.stride * iSubStep;
  const ptrdiff_t iStrideSrc2 = rcDtParam.cur.stride * iSubStep;

  uint32_t uiSum = 0;

#  ifdef USE_AVX2
  if (vext >= AVX2) {
    // Do for width that multiple of 16
    __m256i vone = _mm256_set1_epi16(1);
    __m256i vsum32 = _mm256_setzero_si256();

    for (int iY = 0; iY < iRows; iY += iSubStep) {
      __m256i vsrc1 = _mm256_loadu_si256((__m256i*)(pSrc1));
      __m256i vsrc2 = _mm256_loadu_si256((__m256i*)(pSrc2));
      __m256i vsum16 = _mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2));

      for (int iX = 16; iX < iWidth; iX += 16) {
        vsrc1 = _mm256_loadu_si256((__m256i*)(&pSrc1[iX]));
        vsrc2 = _mm256_loadu_si256((__m256i*)(&pSrc2[iX]));
        vsum16 = _mm256_add_epi16(vsum16, _mm256_abs_epi16(_mm256_sub_epi16(vsrc1, vsrc2)));
      }

      __m256i vsumtemp = _mm256_madd_epi16(vsum16, vone);
      vsum32 = _mm256_add_epi32(vsum32, vsumtemp);

      pSrc1 += iStrideSrc1;
      pSrc2 += iStrideSrc2;
    }

    vsum32 = _mm256_hadd_epi32(vsum32, vone);
    vsum32 = _mm256_hadd_epi32(vsum32, vone);
    uiSum = _mm_cvtsi128_si32(_mm256_castsi256_si128(vsum32)) + _mm_cvtsi128_si32(_mm256_extracti128_si256(vsum32, 1));
  } else
#  endif
  {
    // For width that multiple of 16
    __m128i vone = _mm_set1_epi16(1);
    __m128i vsum32 = _mm_setzero_si128();

    for (int iY = 0; iY < iRows; iY += iSubStep) {
      __m128i vsrc1 = _mm_loadu_si128((const __m128i*)(pSrc1));
      __m128i vsrc2 = _mm_loadu_si128((const __m128i*)(pSrc2));
      __m128i vsum16 = _mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2));

      vsrc1 = _mm_loadu_si128((const __m128i*)(&pSrc1[8]));
      vsrc2 = _mm_loadu_si128((const __m128i*)(&pSrc2[8]));
      vsum16 = _mm_add_epi16(vsum16, _mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2)));

      for (int iX = 16; iX < iWidth; iX += 16) {
        vsrc1 = _mm_loadu_si128((const __m128i*)(&pSrc1[iX]));
        vsrc2 = _mm_loadu_si128((const __m128i*)(&pSrc2[iX]));
        vsum16 = _mm_add_epi16(vsum16, _mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2)));

        vsrc1 = _mm_loadu_si128((const __m128i*)(&pSrc1[iX + 8]));
        vsrc2 = _mm_loadu_si128((const __m128i*)(&pSrc2[iX + 8]));
        vsum16 = _mm_add_epi16(vsum16, _mm_abs_epi16(_mm_sub_epi16(vsrc1, vsrc2)));
      }

      __m128i vsumtemp = _mm_madd_epi16(vsum16, vone);
      vsum32 = _mm_add_epi32(vsum32, vsumtemp);

      pSrc1 += iStrideSrc1;
      pSrc2 += iStrideSrc2;
    }
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    vsum32 = _mm_hadd_epi32(vsum32, vone);
    uiSum = _mm_cvtsi128_si32(vsum32);
  }

  uiSum <<= iSubShift;
  return uiSum >> DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth);
}

template <X86_VEXT vext, bool isCalCentrePos>
void xGetSADX5_4xN_SIMDImp(const DistParam& rcDtParam, Distortion* cost) {
  int i;
  const Pel* piOrg = rcDtParam.org.buf;
  const Pel* piCur = rcDtParam.cur.buf - 4;
  int height = rcDtParam.org.height;
  int iSubShift = rcDtParam.subShift;
  int iSubStep = (1 << iSubShift);
  ptrdiff_t iStrideCur = rcDtParam.cur.stride * iSubStep;
  ptrdiff_t iStrideOrg = rcDtParam.org.stride * iSubStep;

  __m128i sum0 = _mm_setzero_si128();
  __m128i sum1 = _mm_setzero_si128();
  __m128i sum2 = _mm_setzero_si128();
  __m128i sum3 = _mm_setzero_si128();
  __m128i sum4 = _mm_setzero_si128();

  for (i = 0; i < height; i += iSubStep) {
    __m128i s0 = _mm_loadu_si128((__m128i*)piOrg);
    __m128i s1 = _mm_loadu_si128((__m128i*)piCur);

    __m128i org0, org1, org2, org3, org4;
    org0 = s0;
    org1 = _mm_srli_si128(s0, 1);
    if (isCalCentrePos) org2 = _mm_srli_si128(s0, 2);
    org3 = _mm_srli_si128(s0, 3);
    org4 = _mm_srli_si128(s0, 4);

    __m128i cur0, cur1, cur2, cur3, cur4;
    cur4 = s1;
    cur0 = _mm_srli_si128(s1, 4);
    cur1 = _mm_srli_si128(s1, 3);
    if (isCalCentrePos) cur2 = _mm_srli_si128(s1, 2);
    cur3 = _mm_srli_si128(s1, 1);

    __m128i diff0, diff1, diff2, diff3, diff4;
    diff0 = _mm_sub_epi16(org0, cur0);
    diff1 = _mm_sub_epi16(org1, cur1);
    if (isCalCentrePos) diff2 = _mm_sub_epi16(org2, cur2);
    diff3 = _mm_sub_epi16(org3, cur3);
    diff4 = _mm_sub_epi16(org4, cur4);

    diff0 = _mm_abs_epi16(diff0);
    diff1 = _mm_abs_epi16(diff1);
    if (isCalCentrePos) diff2 = _mm_abs_epi16(diff2);
    diff3 = _mm_abs_epi16(diff3);
    diff4 = _mm_abs_epi16(diff4);

    diff0 = _mm_cvtepi16_epi32(diff0);
    diff1 = _mm_cvtepi16_epi32(diff1);
    if (isCalCentrePos) diff2 = _mm_cvtepi16_epi32(diff2);
    diff3 = _mm_cvtepi16_epi32(diff3);
    diff4 = _mm_cvtepi16_epi32(diff4);

    sum0 = _mm_add_epi32(sum0, diff0);
    sum1 = _mm_add_epi32(sum1, diff1);
    if (isCalCentrePos) sum2 = _mm_add_epi32(sum2, diff2);
    sum3 = _mm_add_epi32(sum3, diff3);
    sum4 = _mm_add_epi32(sum4, diff4);

    INCY(piOrg, iStrideOrg);
    INCY(piCur, iStrideCur);
  }

  sum0 = _mm_hadd_epi32(sum0, sum1);
  sum3 = _mm_hadd_epi32(sum3, sum4);
  if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

  sum0 = _mm_hadd_epi32(sum0, sum3);
  if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

  sum0 = _mm_slli_epi32(sum0, iSubShift);
  if (isCalCentrePos) sum2 = _mm_slli_epi32(sum2, iSubShift);

  sum0 = _mm_srli_epi32(sum0, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));
  if (isCalCentrePos) sum2 = _mm_srli_epi32(sum2, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));

  cost[0] = (_mm_cvtsi128_si32(sum0));
  cost[1] = (_mm_extract_epi32(sum0, 1));
  if (isCalCentrePos) cost[2] = (_mm_cvtsi128_si32(sum2));
  cost[3] = (_mm_extract_epi32(sum0, 2));
  cost[4] = (_mm_extract_epi32(sum0, 3));
}

template <X86_VEXT vext>
void RdCost::xGetSADX5_4xN_SIMD(const DistParam& rcDtParam, Distortion* cost, bool isCalCentrePos) {
  if (isCalCentrePos)
    xGetSADX5_4xN_SIMDImp<vext, true>(rcDtParam, cost);
  else
    xGetSADX5_4xN_SIMDImp<vext, false>(rcDtParam, cost);
}

template <X86_VEXT vext, bool isCalCentrePos>
void xGetSADX5_8xN_SIMDImp(const DistParam& rcDtParam, Distortion* cost) {
  int i;
  const Pel* piOrg = rcDtParam.org.buf;
  const Pel* piCur = rcDtParam.cur.buf - 4;
  int height = rcDtParam.org.height;
  int iSubShift = rcDtParam.subShift;
  int iSubStep = (1 << iSubShift);
  ptrdiff_t iStrideCur = rcDtParam.cur.stride * iSubStep;
  ptrdiff_t iStrideOrg = rcDtParam.org.stride * iSubStep;

  __m128i sum0 = _mm_setzero_si128();
  __m128i sum1 = _mm_setzero_si128();
  __m128i sum2 = _mm_setzero_si128();
  __m128i sum3 = _mm_setzero_si128();
  __m128i sum4 = _mm_setzero_si128();

  __m128i vone = _mm_set1_epi16(1);
  for (i = 0; i < height; i += iSubStep) {
    __m128i s0 = _mm_loadu_si128((__m128i*)piOrg);
    __m128i s1 = _mm_loadu_si128((__m128i*)piCur);
    __m128i s2 = _mm_loadl_epi64((__m128i*)(piOrg + 8));
    __m128i s3 = _mm_loadl_epi64((__m128i*)(piCur + 8));

    __m128i org0, org1, org2, org3, org4;
    org0 = s0;
    org1 = _mm_alignr_epi8(s2, s0, 2);
    if (isCalCentrePos) org2 = _mm_alignr_epi8(s2, s0, 4);
    org3 = _mm_alignr_epi8(s2, s0, 6);
    org4 = _mm_alignr_epi8(s2, s0, 8);

    __m128i cur0, cur1, cur2, cur3, cur4;
    cur4 = s1;
    cur0 = _mm_alignr_epi8(s3, s1, 8);
    cur1 = _mm_alignr_epi8(s3, s1, 6);
    if (isCalCentrePos) cur2 = _mm_alignr_epi8(s3, s1, 4);
    cur3 = _mm_alignr_epi8(s3, s1, 2);

    __m128i diff0, diff1, diff2, diff3, diff4;
    diff0 = _mm_sub_epi16(org0, cur0);
    diff1 = _mm_sub_epi16(org1, cur1);
    if (isCalCentrePos) diff2 = _mm_sub_epi16(org2, cur2);
    diff3 = _mm_sub_epi16(org3, cur3);
    diff4 = _mm_sub_epi16(org4, cur4);

    diff0 = _mm_abs_epi16(diff0);
    diff1 = _mm_abs_epi16(diff1);
    if (isCalCentrePos) diff2 = _mm_abs_epi16(diff2);
    diff3 = _mm_abs_epi16(diff3);
    diff4 = _mm_abs_epi16(diff4);

    diff0 = _mm_madd_epi16(diff0, vone);
    diff1 = _mm_madd_epi16(diff1, vone);
    if (isCalCentrePos) diff2 = _mm_madd_epi16(diff2, vone);
    diff3 = _mm_madd_epi16(diff3, vone);
    diff4 = _mm_madd_epi16(diff4, vone);

    sum0 = _mm_add_epi32(sum0, diff0);
    sum1 = _mm_add_epi32(sum1, diff1);
    if (isCalCentrePos) sum2 = _mm_add_epi32(sum2, diff2);
    sum3 = _mm_add_epi32(sum3, diff3);
    sum4 = _mm_add_epi32(sum4, diff4);

    INCY(piOrg, iStrideOrg);
    INCY(piCur, iStrideCur);
  }

  sum0 = _mm_hadd_epi32(sum0, sum1);
  sum3 = _mm_hadd_epi32(sum3, sum4);
  if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

  sum0 = _mm_hadd_epi32(sum0, sum3);
  if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

  sum0 = _mm_slli_epi32(sum0, iSubShift);
  if (isCalCentrePos) sum2 = _mm_slli_epi32(sum2, iSubShift);

  sum0 = _mm_srli_epi32(sum0, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));
  if (isCalCentrePos) sum2 = _mm_srli_epi32(sum2, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));

  cost[0] = (_mm_cvtsi128_si32(sum0));
  cost[1] = (_mm_extract_epi32(sum0, 1));
  if (isCalCentrePos) cost[2] = (_mm_cvtsi128_si32(sum2));
  cost[3] = (_mm_extract_epi32(sum0, 2));
  cost[4] = (_mm_extract_epi32(sum0, 3));
}

template <X86_VEXT vext>
void RdCost::xGetSADX5_8xN_SIMD(const DistParam& rcDtParam, Distortion* cost, bool isCalCentrePos) {
  if (isCalCentrePos)
    xGetSADX5_8xN_SIMDImp<vext, true>(rcDtParam, cost);
  else
    xGetSADX5_8xN_SIMDImp<vext, false>(rcDtParam, cost);
}

template <int iWidth, X86_VEXT vext, bool isCalCentrePos>
void xGetSADX5_16NxN_SIMDImp(const DistParam& rcDtParam, Distortion* cost) {
  int i, j;
  const Pel* piOrg = rcDtParam.org.buf;
  const Pel* piCur = rcDtParam.cur.buf - 4;
  int height = rcDtParam.org.height;
  int width = rcDtParam.org.width;
  int iSubShift = rcDtParam.subShift;
  int iSubStep = (1 << iSubShift);
  ptrdiff_t iStrideCur = rcDtParam.cur.stride * iSubStep;
  ptrdiff_t iStrideOrg = rcDtParam.org.stride * iSubStep;

#  ifdef USE_AVX2
  if (vext >= AVX2) {
    __m256i sum0 = _mm256_setzero_si256();
    __m256i sum1 = _mm256_setzero_si256();
    __m256i sum2 = _mm256_setzero_si256();
    __m256i sum3 = _mm256_setzero_si256();
    __m256i sum4 = _mm256_setzero_si256();

    __m256i vone = _mm256_set1_epi16(1);
    for (i = 0; i < height; i += iSubStep) {
      __m256i sumTmp0 = _mm256_setzero_si256();
      __m256i sumTmp1 = _mm256_setzero_si256();
      __m256i sumTmp2 = _mm256_setzero_si256();
      __m256i sumTmp3 = _mm256_setzero_si256();
      __m256i sumTmp4 = _mm256_setzero_si256();
      for (j = 0; j < width; j += 16) {
        __m256i s0 = _mm256_loadu_si256((__m256i*)piOrg);
        __m256i s1 = _mm256_loadu_si256((__m256i*)piCur);
        __m256i s2 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piOrg + 16)));
        __m256i s3 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piCur + 16)));
        s2 = _mm256_permute2x128_si256(s0, s2, 0x21);
        s3 = _mm256_permute2x128_si256(s1, s3, 0x21);

        __m256i org0, org1, org2, org3, org4;
        org0 = s0;
        org1 = _mm256_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm256_alignr_epi8(s2, s0, 4);
        org3 = _mm256_alignr_epi8(s2, s0, 6);
        org4 = _mm256_alignr_epi8(s2, s0, 8);

        __m256i cur0, cur1, cur2, cur3, cur4;
        cur4 = s1;
        cur0 = _mm256_alignr_epi8(s3, s1, 8);
        cur1 = _mm256_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm256_alignr_epi8(s3, s1, 4);
        cur3 = _mm256_alignr_epi8(s3, s1, 2);

        __m256i diff0, diff1, diff2, diff3, diff4;
        diff0 = _mm256_sub_epi16(org0, cur0);
        diff1 = _mm256_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm256_sub_epi16(org2, cur2);
        diff3 = _mm256_sub_epi16(org3, cur3);
        diff4 = _mm256_sub_epi16(org4, cur4);

        diff0 = _mm256_abs_epi16(diff0);
        diff1 = _mm256_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm256_abs_epi16(diff2);
        diff3 = _mm256_abs_epi16(diff3);
        diff4 = _mm256_abs_epi16(diff4);

        sumTmp0 = _mm256_add_epi16(sumTmp0, diff0);
        sumTmp1 = _mm256_add_epi16(sumTmp1, diff1);
        if (isCalCentrePos) sumTmp2 = _mm256_add_epi16(sumTmp2, diff2);
        sumTmp3 = _mm256_add_epi16(sumTmp3, diff3);
        sumTmp4 = _mm256_add_epi16(sumTmp4, diff4);
      }

      sumTmp0 = _mm256_madd_epi16(sumTmp0, vone);
      sumTmp1 = _mm256_madd_epi16(sumTmp1, vone);
      if (isCalCentrePos) sumTmp2 = _mm256_madd_epi16(sumTmp2, vone);
      sumTmp3 = _mm256_madd_epi16(sumTmp3, vone);
      sumTmp4 = _mm256_madd_epi16(sumTmp4, vone);

      sum0 = _mm256_add_epi32(sum0, sumTmp0);
      sum1 = _mm256_add_epi32(sum1, sumTmp1);
      if (isCalCentrePos) sum2 = _mm256_add_epi32(sum2, sumTmp2);
      sum3 = _mm256_add_epi32(sum3, sumTmp3);
      sum4 = _mm256_add_epi32(sum4, sumTmp4);

      INCY(piOrg, iStrideOrg);
      INCY(piCur, iStrideCur);
    }

    sum0 = _mm256_hadd_epi32(sum0, sum1);
    sum3 = _mm256_hadd_epi32(sum3, sum4);
    if (isCalCentrePos) sum2 = _mm256_hadd_epi32(sum2, sum2);

    sum0 = _mm256_hadd_epi32(sum0, sum3);
    if (isCalCentrePos) sum2 = _mm256_hadd_epi32(sum2, sum2);

    __m128i sum0134 = _mm_add_epi32(_mm256_castsi256_si128(sum0), _mm256_extracti128_si256(sum0, 1));

    sum0134 = _mm_slli_epi32(sum0134, iSubShift);

    sum0134 = _mm_srli_epi32(sum0134, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));

    cost[0] = (_mm_cvtsi128_si32(sum0134));
    cost[1] = (_mm_extract_epi32(sum0134, 1));
    if (isCalCentrePos) {
      int tmp = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum2)) + _mm256_extract_epi32(sum2, 4);
      tmp <<= iSubShift;
      tmp >>= (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth)));
      cost[2] = tmp;
    }
    cost[3] = (_mm_extract_epi32(sum0134, 2));
    cost[4] = (_mm_extract_epi32(sum0134, 3));
  } else
#  endif
  {
    __m128i sum0 = _mm_setzero_si128();
    __m128i sum1 = _mm_setzero_si128();
    __m128i sum2 = _mm_setzero_si128();
    __m128i sum3 = _mm_setzero_si128();
    __m128i sum4 = _mm_setzero_si128();

    __m128i vone = _mm_set1_epi16(1);
    for (i = 0; i < height; i += iSubStep) {
      __m128i sumTmp0 = _mm_setzero_si128();
      __m128i sumTmp1 = _mm_setzero_si128();
      __m128i sumTmp2 = _mm_setzero_si128();
      __m128i sumTmp3 = _mm_setzero_si128();
      __m128i sumTmp4 = _mm_setzero_si128();
      for (j = 0; j < width; j += 8) {
        __m128i s0 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(piOrg + 0));
        __m128i s1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(piCur + 0));
        __m128i s2 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(piOrg + 8));
        __m128i s3 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(piCur + 8));

        __m128i org0, org1, org2, org3, org4;
        org0 = s0;
        org1 = _mm_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm_alignr_epi8(s2, s0, 4);
        org3 = _mm_alignr_epi8(s2, s0, 6);
        org4 = _mm_alignr_epi8(s2, s0, 8);

        __m128i cur0, cur1, cur2, cur3, cur4;
        cur4 = s1;
        cur0 = _mm_alignr_epi8(s3, s1, 8);
        cur1 = _mm_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm_alignr_epi8(s3, s1, 4);
        cur3 = _mm_alignr_epi8(s3, s1, 2);

        __m128i diff0, diff1, diff2, diff3, diff4;
        diff0 = _mm_sub_epi16(org0, cur0);
        diff1 = _mm_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm_sub_epi16(org2, cur2);
        diff3 = _mm_sub_epi16(org3, cur3);
        diff4 = _mm_sub_epi16(org4, cur4);

        diff0 = _mm_abs_epi16(diff0);
        diff1 = _mm_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm_abs_epi16(diff2);
        diff3 = _mm_abs_epi16(diff3);
        diff4 = _mm_abs_epi16(diff4);

        sumTmp0 = _mm_add_epi16(sumTmp0, diff0);
        sumTmp1 = _mm_add_epi16(sumTmp1, diff1);
        if (isCalCentrePos) sumTmp2 = _mm_add_epi16(sumTmp2, diff2);
        sumTmp3 = _mm_add_epi16(sumTmp3, diff3);
        sumTmp4 = _mm_add_epi16(sumTmp4, diff4);
      }

      sumTmp0 = _mm_madd_epi16(sumTmp0, vone);
      sumTmp1 = _mm_madd_epi16(sumTmp1, vone);
      if (isCalCentrePos) sumTmp2 = _mm_madd_epi16(sumTmp2, vone);
      sumTmp3 = _mm_madd_epi16(sumTmp3, vone);
      sumTmp4 = _mm_madd_epi16(sumTmp4, vone);

      sum0 = _mm_add_epi32(sum0, sumTmp0);
      sum1 = _mm_add_epi32(sum1, sumTmp1);
      if (isCalCentrePos) sum2 = _mm_add_epi32(sum2, sumTmp2);
      sum3 = _mm_add_epi32(sum3, sumTmp3);
      sum4 = _mm_add_epi32(sum4, sumTmp4);

      INCY(piOrg, iStrideOrg);
      INCY(piCur, iStrideCur);
    }

    sum0 = _mm_hadd_epi32(sum0, sum1);
    sum3 = _mm_hadd_epi32(sum3, sum4);
    if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

    sum0 = _mm_hadd_epi32(sum0, sum3);
    if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

    sum0 = _mm_slli_epi32(sum0, iSubShift);
    if (isCalCentrePos) sum2 = _mm_slli_epi32(sum2, iSubShift);

    sum0 = _mm_srli_epi32(sum0, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));
    if (isCalCentrePos) sum2 = _mm_srli_epi32(sum2, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));

    cost[0] = (_mm_cvtsi128_si32(sum0));
    cost[1] = (_mm_extract_epi32(sum0, 1));
    if (isCalCentrePos) cost[2] = (_mm_cvtsi128_si32(sum2));
    cost[3] = (_mm_extract_epi32(sum0, 2));
    cost[4] = (_mm_extract_epi32(sum0, 3));
  }
}

template <int iWidth, X86_VEXT vext>
void RdCost::xGetSADX5_16NxN_SIMD(const DistParam& rcDtParam, Distortion* cost, bool isCalCentrePos) {
  if (isCalCentrePos)
    xGetSADX5_16NxN_SIMDImp<iWidth, vext, true>(rcDtParam, cost);
  else
    xGetSADX5_16NxN_SIMDImp<iWidth, vext, false>(rcDtParam, cost);
}

template <X86_VEXT vext, bool isCalCentrePos>
void xGetSADX5_16xN_SIMDImp(const DistParam& rcDtParam, Distortion* cost) {
  int i, j;
  const Pel* piOrg = rcDtParam.org.buf;
  const Pel* piCur = rcDtParam.cur.buf - 4;
  int height = rcDtParam.org.height;
  int iSubShift = rcDtParam.subShift;
  int iSubStep = (1 << iSubShift);
  ptrdiff_t iStrideCur = rcDtParam.cur.stride * iSubStep;
  ptrdiff_t iStrideOrg = rcDtParam.org.stride * iSubStep;

#  ifdef USE_AVX2
  if (vext >= AVX2) {
    __m256i sum0 = _mm256_setzero_si256();
    __m256i sum1 = _mm256_setzero_si256();
    __m256i sum2 = _mm256_setzero_si256();
    __m256i sum3 = _mm256_setzero_si256();
    __m256i sum4 = _mm256_setzero_si256();

    __m256i vone = _mm256_set1_epi16(1);

    if (height == 16 && iSubShift == 1) {
      __m256i sumTmp0, sumTmp1, sumTmp2, sumTmp3, sumTmp4;
      for (int i = 0; i < 2; i++) {
        __m256i s0 = _mm256_loadu_si256((__m256i*)piOrg);
        __m256i s1 = _mm256_loadu_si256((__m256i*)piCur);
        __m256i s2 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piOrg + 16)));
        __m256i s3 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piCur + 16)));
        s2 = _mm256_permute2x128_si256(s0, s2, 0x21);
        s3 = _mm256_permute2x128_si256(s1, s3, 0x21);

        INCY(piOrg, iStrideOrg);
        INCY(piCur, iStrideCur);

        __m256i org0, org1, org2, org3, org4;
        org0 = s0;
        org1 = _mm256_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm256_alignr_epi8(s2, s0, 4);
        org3 = _mm256_alignr_epi8(s2, s0, 6);
        org4 = _mm256_alignr_epi8(s2, s0, 8);

        __m256i cur0, cur1, cur2, cur3, cur4;
        cur4 = s1;
        cur0 = _mm256_alignr_epi8(s3, s1, 8);
        cur1 = _mm256_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm256_alignr_epi8(s3, s1, 4);
        cur3 = _mm256_alignr_epi8(s3, s1, 2);

        __m256i diff0, diff1, diff2, diff3, diff4;
        diff0 = _mm256_sub_epi16(org0, cur0);
        diff1 = _mm256_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm256_sub_epi16(org2, cur2);
        diff3 = _mm256_sub_epi16(org3, cur3);
        diff4 = _mm256_sub_epi16(org4, cur4);

        sumTmp0 = _mm256_abs_epi16(diff0);
        sumTmp1 = _mm256_abs_epi16(diff1);
        if (isCalCentrePos) sumTmp2 = _mm256_abs_epi16(diff2);
        sumTmp3 = _mm256_abs_epi16(diff3);
        sumTmp4 = _mm256_abs_epi16(diff4);

        s0 = _mm256_loadu_si256((__m256i*)piOrg);
        s1 = _mm256_loadu_si256((__m256i*)piCur);
        s2 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piOrg + 16)));
        s3 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piCur + 16)));
        s2 = _mm256_permute2x128_si256(s0, s2, 0x21);
        s3 = _mm256_permute2x128_si256(s1, s3, 0x21);

        INCY(piOrg, iStrideOrg);
        INCY(piCur, iStrideCur);

        org0 = s0;
        org1 = _mm256_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm256_alignr_epi8(s2, s0, 4);
        org3 = _mm256_alignr_epi8(s2, s0, 6);
        org4 = _mm256_alignr_epi8(s2, s0, 8);

        cur4 = s1;
        cur0 = _mm256_alignr_epi8(s3, s1, 8);
        cur1 = _mm256_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm256_alignr_epi8(s3, s1, 4);
        cur3 = _mm256_alignr_epi8(s3, s1, 2);

        diff0 = _mm256_sub_epi16(org0, cur0);
        diff1 = _mm256_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm256_sub_epi16(org2, cur2);
        diff3 = _mm256_sub_epi16(org3, cur3);
        diff4 = _mm256_sub_epi16(org4, cur4);

        diff0 = _mm256_abs_epi16(diff0);
        diff1 = _mm256_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm256_abs_epi16(diff2);
        diff3 = _mm256_abs_epi16(diff3);
        diff4 = _mm256_abs_epi16(diff4);

        sumTmp0 = _mm256_add_epi16(diff0, sumTmp0);
        sumTmp1 = _mm256_add_epi16(diff1, sumTmp1);
        if (isCalCentrePos) sumTmp2 = _mm256_add_epi16(diff2, sumTmp2);
        sumTmp3 = _mm256_add_epi16(diff3, sumTmp3);
        sumTmp4 = _mm256_add_epi16(diff4, sumTmp4);

        s0 = _mm256_loadu_si256((__m256i*)piOrg);
        s1 = _mm256_loadu_si256((__m256i*)piCur);
        s2 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piOrg + 16)));
        s3 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piCur + 16)));
        s2 = _mm256_permute2x128_si256(s0, s2, 0x21);
        s3 = _mm256_permute2x128_si256(s1, s3, 0x21);

        INCY(piOrg, iStrideOrg);
        INCY(piCur, iStrideCur);

        org0 = s0;
        org1 = _mm256_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm256_alignr_epi8(s2, s0, 4);
        org3 = _mm256_alignr_epi8(s2, s0, 6);
        org4 = _mm256_alignr_epi8(s2, s0, 8);

        cur4 = s1;
        cur0 = _mm256_alignr_epi8(s3, s1, 8);
        cur1 = _mm256_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm256_alignr_epi8(s3, s1, 4);
        cur3 = _mm256_alignr_epi8(s3, s1, 2);

        diff0 = _mm256_sub_epi16(org0, cur0);
        diff1 = _mm256_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm256_sub_epi16(org2, cur2);
        diff3 = _mm256_sub_epi16(org3, cur3);
        diff4 = _mm256_sub_epi16(org4, cur4);

        diff0 = _mm256_abs_epi16(diff0);
        diff1 = _mm256_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm256_abs_epi16(diff2);
        diff3 = _mm256_abs_epi16(diff3);
        diff4 = _mm256_abs_epi16(diff4);

        sumTmp0 = _mm256_add_epi16(diff0, sumTmp0);
        sumTmp1 = _mm256_add_epi16(diff1, sumTmp1);
        if (isCalCentrePos) sumTmp2 = _mm256_add_epi16(diff2, sumTmp2);
        sumTmp3 = _mm256_add_epi16(diff3, sumTmp3);
        sumTmp4 = _mm256_add_epi16(diff4, sumTmp4);

        s0 = _mm256_loadu_si256((__m256i*)piOrg);
        s1 = _mm256_loadu_si256((__m256i*)piCur);
        s2 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piOrg + 16)));
        s3 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piCur + 16)));
        s2 = _mm256_permute2x128_si256(s0, s2, 0x21);
        s3 = _mm256_permute2x128_si256(s1, s3, 0x21);

        INCY(piOrg, iStrideOrg);
        INCY(piCur, iStrideCur);

        org0 = s0;
        org1 = _mm256_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm256_alignr_epi8(s2, s0, 4);
        org3 = _mm256_alignr_epi8(s2, s0, 6);
        org4 = _mm256_alignr_epi8(s2, s0, 8);

        cur4 = s1;
        cur0 = _mm256_alignr_epi8(s3, s1, 8);
        cur1 = _mm256_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm256_alignr_epi8(s3, s1, 4);
        cur3 = _mm256_alignr_epi8(s3, s1, 2);

        diff0 = _mm256_sub_epi16(org0, cur0);
        diff1 = _mm256_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm256_sub_epi16(org2, cur2);
        diff3 = _mm256_sub_epi16(org3, cur3);
        diff4 = _mm256_sub_epi16(org4, cur4);

        diff0 = _mm256_abs_epi16(diff0);
        diff1 = _mm256_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm256_abs_epi16(diff2);
        diff3 = _mm256_abs_epi16(diff3);
        diff4 = _mm256_abs_epi16(diff4);

        sumTmp0 = _mm256_add_epi16(diff0, sumTmp0);
        sumTmp1 = _mm256_add_epi16(diff1, sumTmp1);
        if (isCalCentrePos) sumTmp2 = _mm256_add_epi16(diff2, sumTmp2);
        sumTmp3 = _mm256_add_epi16(diff3, sumTmp3);
        sumTmp4 = _mm256_add_epi16(diff4, sumTmp4);

        sumTmp0 = _mm256_madd_epi16(sumTmp0, vone);
        sumTmp1 = _mm256_madd_epi16(sumTmp1, vone);
        if (isCalCentrePos) sumTmp2 = _mm256_madd_epi16(sumTmp2, vone);
        sumTmp3 = _mm256_madd_epi16(sumTmp3, vone);
        sumTmp4 = _mm256_madd_epi16(sumTmp4, vone);

        sum0 = _mm256_add_epi32(sum0, sumTmp0);
        sum1 = _mm256_add_epi32(sum1, sumTmp1);
        if (isCalCentrePos) sum2 = _mm256_add_epi32(sum2, sumTmp2);
        sum3 = _mm256_add_epi32(sum3, sumTmp3);
        sum4 = _mm256_add_epi32(sum4, sumTmp4);
      }
    } else {
      for (i = 0; i < height; i += iSubStep) {
        __m256i s0 = _mm256_loadu_si256((__m256i*)piOrg);
        __m256i s1 = _mm256_loadu_si256((__m256i*)piCur);
        __m256i s2 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piOrg + 16)));
        __m256i s3 = _mm256_castsi128_si256(_mm_loadl_epi64((__m128i*)(piCur + 16)));
        s2 = _mm256_permute2x128_si256(s0, s2, 0x21);
        s3 = _mm256_permute2x128_si256(s1, s3, 0x21);

        __m256i org0, org1, org2, org3, org4;
        org0 = s0;
        org1 = _mm256_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm256_alignr_epi8(s2, s0, 4);
        org3 = _mm256_alignr_epi8(s2, s0, 6);
        org4 = _mm256_alignr_epi8(s2, s0, 8);

        __m256i cur0, cur1, cur2, cur3, cur4;
        cur4 = s1;
        cur0 = _mm256_alignr_epi8(s3, s1, 8);
        cur1 = _mm256_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm256_alignr_epi8(s3, s1, 4);
        cur3 = _mm256_alignr_epi8(s3, s1, 2);

        __m256i diff0, diff1, diff2, diff3, diff4;
        diff0 = _mm256_sub_epi16(org0, cur0);
        diff1 = _mm256_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm256_sub_epi16(org2, cur2);
        diff3 = _mm256_sub_epi16(org3, cur3);
        diff4 = _mm256_sub_epi16(org4, cur4);

        diff0 = _mm256_abs_epi16(diff0);
        diff1 = _mm256_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm256_abs_epi16(diff2);
        diff3 = _mm256_abs_epi16(diff3);
        diff4 = _mm256_abs_epi16(diff4);

        diff0 = _mm256_madd_epi16(diff0, vone);
        diff1 = _mm256_madd_epi16(diff1, vone);
        if (isCalCentrePos) diff2 = _mm256_madd_epi16(diff2, vone);
        diff3 = _mm256_madd_epi16(diff3, vone);
        diff4 = _mm256_madd_epi16(diff4, vone);

        sum0 = _mm256_add_epi32(sum0, diff0);
        sum1 = _mm256_add_epi32(sum1, diff1);
        if (isCalCentrePos) sum2 = _mm256_add_epi32(sum2, diff2);
        sum3 = _mm256_add_epi32(sum3, diff3);
        sum4 = _mm256_add_epi32(sum4, diff4);

        INCY(piOrg, iStrideOrg);
        INCY(piCur, iStrideCur);
      }
    }

    sum0 = _mm256_hadd_epi32(sum0, sum1);
    sum3 = _mm256_hadd_epi32(sum3, sum4);
    if (isCalCentrePos) sum2 = _mm256_hadd_epi32(sum2, sum2);

    sum0 = _mm256_hadd_epi32(sum0, sum3);
    if (isCalCentrePos) sum2 = _mm256_hadd_epi32(sum2, sum2);

    __m128i sum0134 = _mm_add_epi32(_mm256_castsi256_si128(sum0), _mm256_extracti128_si256(sum0, 1));

    sum0134 = _mm_slli_epi32(sum0134, iSubShift);

    sum0134 = _mm_srli_epi32(sum0134, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));

    cost[0] = (_mm_cvtsi128_si32(sum0134));
    cost[1] = (_mm_extract_epi32(sum0134, 1));
    if (isCalCentrePos) {
      int tmp = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum2)) + _mm256_extract_epi32(sum2, 4);
      tmp <<= iSubShift;
      tmp >>= (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth)));
      cost[2] = tmp;
    }
    cost[3] = (_mm_extract_epi32(sum0134, 2));
    cost[4] = (_mm_extract_epi32(sum0134, 3));
  } else
#  endif
  {
    __m128i sum0 = _mm_setzero_si128();
    __m128i sum1 = _mm_setzero_si128();
    __m128i sum2 = _mm_setzero_si128();
    __m128i sum3 = _mm_setzero_si128();
    __m128i sum4 = _mm_setzero_si128();

    __m128i vone = _mm_set1_epi16(1);
    for (i = 0; i < height; i += iSubStep) {
      __m128i sumTmp0 = _mm_setzero_si128();
      __m128i sumTmp1 = _mm_setzero_si128();
      __m128i sumTmp2 = _mm_setzero_si128();
      __m128i sumTmp3 = _mm_setzero_si128();
      __m128i sumTmp4 = _mm_setzero_si128();
      for (j = 0; j < 16; j += 8) {
        __m128i s0 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(piOrg + j + 0));
        __m128i s1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(piCur + j + 0));
        __m128i s2 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(piOrg + j + 8));
        __m128i s3 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(piCur + j + 8));

        __m128i org0, org1, org2, org3, org4;
        org0 = s0;
        org1 = _mm_alignr_epi8(s2, s0, 2);
        if (isCalCentrePos) org2 = _mm_alignr_epi8(s2, s0, 4);
        org3 = _mm_alignr_epi8(s2, s0, 6);
        org4 = _mm_alignr_epi8(s2, s0, 8);

        __m128i cur0, cur1, cur2, cur3, cur4;
        cur4 = s1;
        cur0 = _mm_alignr_epi8(s3, s1, 8);
        cur1 = _mm_alignr_epi8(s3, s1, 6);
        if (isCalCentrePos) cur2 = _mm_alignr_epi8(s3, s1, 4);
        cur3 = _mm_alignr_epi8(s3, s1, 2);

        __m128i diff0, diff1, diff2, diff3, diff4;
        diff0 = _mm_sub_epi16(org0, cur0);
        diff1 = _mm_sub_epi16(org1, cur1);
        if (isCalCentrePos) diff2 = _mm_sub_epi16(org2, cur2);
        diff3 = _mm_sub_epi16(org3, cur3);
        diff4 = _mm_sub_epi16(org4, cur4);

        diff0 = _mm_abs_epi16(diff0);
        diff1 = _mm_abs_epi16(diff1);
        if (isCalCentrePos) diff2 = _mm_abs_epi16(diff2);
        diff3 = _mm_abs_epi16(diff3);
        diff4 = _mm_abs_epi16(diff4);

        sumTmp0 = _mm_add_epi16(sumTmp0, diff0);
        sumTmp1 = _mm_add_epi16(sumTmp1, diff1);
        if (isCalCentrePos) sumTmp2 = _mm_add_epi16(sumTmp2, diff2);
        sumTmp3 = _mm_add_epi16(sumTmp3, diff3);
        sumTmp4 = _mm_add_epi16(sumTmp4, diff4);
      }

      sumTmp0 = _mm_madd_epi16(sumTmp0, vone);
      sumTmp1 = _mm_madd_epi16(sumTmp1, vone);
      if (isCalCentrePos) sumTmp2 = _mm_madd_epi16(sumTmp2, vone);
      sumTmp3 = _mm_madd_epi16(sumTmp3, vone);
      sumTmp4 = _mm_madd_epi16(sumTmp4, vone);

      sum0 = _mm_add_epi32(sum0, sumTmp0);
      sum1 = _mm_add_epi32(sum1, sumTmp1);
      if (isCalCentrePos) sum2 = _mm_add_epi32(sum2, sumTmp2);
      sum3 = _mm_add_epi32(sum3, sumTmp3);
      sum4 = _mm_add_epi32(sum4, sumTmp4);

      INCY(piOrg, iStrideOrg);
      INCY(piCur, iStrideCur);
    }

    sum0 = _mm_hadd_epi32(sum0, sum1);
    sum3 = _mm_hadd_epi32(sum3, sum4);
    if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

    sum0 = _mm_hadd_epi32(sum0, sum3);
    if (isCalCentrePos) sum2 = _mm_hadd_epi32(sum2, sum2);

    sum0 = _mm_slli_epi32(sum0, iSubShift);
    if (isCalCentrePos) sum2 = _mm_slli_epi32(sum2, iSubShift);

    sum0 = _mm_srli_epi32(sum0, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));
    if (isCalCentrePos) sum2 = _mm_srli_epi32(sum2, (1 + (DISTORTION_PRECISION_ADJUSTMENT(rcDtParam.bitDepth))));

    cost[0] = (_mm_cvtsi128_si32(sum0));
    cost[1] = (_mm_extract_epi32(sum0, 1));
    if (isCalCentrePos) cost[2] = (_mm_cvtsi128_si32(sum2));
    cost[3] = (_mm_extract_epi32(sum0, 2));
    cost[4] = (_mm_extract_epi32(sum0, 3));
  }
}

template <X86_VEXT vext>
void RdCost::xGetSADX5_16xN_SIMD(const DistParam& rcDtParam, Distortion* cost, bool isCalCentrePos) {
  if (isCalCentrePos)
    xGetSADX5_16xN_SIMDImp<vext, true>(rcDtParam, cost);
  else
    xGetSADX5_16xN_SIMDImp<vext, false>(rcDtParam, cost);
}

template <X86_VEXT vext>
void RdCost::_initRdCostX86(int bitDepth) {
  m_afpDistortFunc[DF_SAD4] = xGetSAD_4xN_SIMD<vext>;
  m_afpDistortFunc[DF_SAD8] = xGetSAD_8xN_SIMD<vext>;
  m_afpDistortFunc[DF_SAD16] = xGetSAD_16xN_SIMD<vext>;
  m_afpDistortFunc[DF_SAD32] = xGetSAD_16NxN_SIMD<32, vext>;
  m_afpDistortFunc[DF_SAD64] = xGetSAD_16NxN_SIMD<64, vext>;
  m_afpDistortFunc[DF_SAD16N] = xGetSAD_16NxN_SIMD<16, vext>;

  m_afpDistortFuncX5[DF_SAD4_X5] = xGetSADX5_4xN_SIMD<vext>;
  m_afpDistortFuncX5[DF_SAD8_X5] = xGetSADX5_8xN_SIMD<vext>;
  m_afpDistortFuncX5[DF_SAD16_X5] = xGetSADX5_16xN_SIMD<vext>;
  m_afpDistortFuncX5[DF_SAD32_X5] = xGetSADX5_16NxN_SIMD<32, vext>;
  m_afpDistortFuncX5[DF_SAD64_X5] = xGetSADX5_16NxN_SIMD<64, vext>;
  m_afpDistortFuncX5[DF_SAD16N_X5] = xGetSADX5_16NxN_SIMD<16, vext>;
}

template void RdCost::_initRdCostX86<SIMDX86>(int bitDepth);

#endif  //#if TARGET_SIMD_X86
//! \}
