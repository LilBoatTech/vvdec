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

/** \file     MatrixIntraPrediction.cpp
\brief    matrix-based intra prediction class
*/

#include "MatrixIntraPrediction.h"
#include "dtrace_next.h"

#include "UnitTools.h"
#include "MipData.h"

namespace Mip {
PredictorMIP::PredictorMIP()
    : m_blockSize(0, 0),
      m_sizeId(0),
      m_reducedBdrySize(0),
      m_reducedPredSize(0),
      m_upsmpFactorHor(0),
      m_upsmpFactorVer(0) {}

void PredictorMIP::deriveBoundaryData(const CPelBuf& src, const Area& block, const int bitDepth) {
  // Step 1: Save block size and calculate dependent values
  initPredBlockParams(block);

  // Step 2: Get the input data (left and top reference samples)
#if ADAPTIVE_BIT_DEPTH
  if (bitDepth <= 8) {
    const Pel8bit* srcPtr = reinterpret_cast<const Pel8bit*>(src.buf);
    ptrdiff_t srcStride = src.stride;
    for (int x = 0; x < block.width; x++) {
      m_refSamplesTop[x] = srcPtr[x + 1];
    }

    //    m_refSamplesLeft.resize(block.height);
    for (int y = 0; y < block.height; y++) {
      m_refSamplesLeft[y] = srcPtr[(y + 1) * srcStride];
    }
  } else {
    for (int x = 0; x < block.width; x++) {
      m_refSamplesTop[x] = src.at(x + 1, 0);
    }

    //    m_refSamplesLeft.resize(block.height);
    for (int y = 0; y < block.height; y++) {
      m_refSamplesLeft[y] = src.at(0, y + 1);
    }
  }
#else
  for (int x = 0; x < block.width; x++) {
    m_refSamplesTop[x] = src.at(x + 1, 0);
  }

  //    m_refSamplesLeft.resize(block.height);
  for (int y = 0; y < block.height; y++) {
    m_refSamplesLeft[y] = src.at(0, y + 1);
  }
#endif

  // Step 3: Compute the reduced boundary via Haar-downsampling (input for the prediction)
  int16_t* const topReduced = m_reducedBoundary;
  boundaryDownsampling1D(topReduced, m_refSamplesTop, block.width, m_reducedBdrySize);

  int16_t* const leftReduced = m_reducedBoundary + m_reducedBdrySize;
  boundaryDownsampling1D(leftReduced, m_refSamplesLeft, block.height, m_reducedBdrySize);

  int16_t* const leftReducedTransposed = m_reducedBoundaryTransposed;
  int16_t* const topReducedTransposed = m_reducedBoundaryTransposed + m_reducedBdrySize;
  for (int x = 0; x < m_reducedBdrySize; x++) {
    topReducedTransposed[x] = topReduced[x];
  }
  for (int y = 0; y < m_reducedBdrySize; y++) {
    leftReducedTransposed[y] = leftReduced[y];
  }

  // Step 4: Rebase the reduced boundary
  const int inputSize = 2 * m_reducedBdrySize;

  m_inputOffset = m_reducedBoundary[0];
  m_inputOffsetTransp = m_reducedBoundaryTransposed[0];

  const bool hasFirstCol = (m_sizeId < 2);
  m_reducedBoundary[0] =
      hasFirstCol ? ((1 << (bitDepth - 1)) - m_inputOffset) : 0;  // first column of matrix not needed for large blocks
  m_reducedBoundaryTransposed[0] = hasFirstCol ? ((1 << (bitDepth - 1)) - m_inputOffsetTransp) : 0;
  for (int i = 1; i < inputSize; i++) {
    m_reducedBoundary[i] -= m_inputOffset;
    m_reducedBoundaryTransposed[i] -= m_inputOffsetTransp;
  }
}

void PredictorMIP::getPrediction(int16_t* const result, const int modeIdx, const bool transpose, const int bitDepth) {
  const bool needUpsampling = (m_upsmpFactorHor > 1) || (m_upsmpFactorVer > 1);

  const uint8_t* matrix = getMatrixData(modeIdx);

  int16_t bufReducedPred[MIP_MAX_REDUCED_OUTPUT_SAMPLES];
  int16_t* const reducedPred = needUpsampling ? bufReducedPred : result;
  const int16_t* const reducedBoundary = transpose ? m_reducedBoundaryTransposed : m_reducedBoundary;
  computeReducedPred(reducedPred, reducedBoundary, matrix, transpose, bitDepth, 2 * m_reducedBdrySize,
                     transpose ? m_inputOffsetTransp : m_inputOffset, m_reducedPredSize, m_sizeId);
  if (needUpsampling) {
    predictionUpsampling(result, reducedPred);
  }
}

void PredictorMIP::initPredBlockParams(const Size& block) {
  m_blockSize = block;
  // init size index
  m_sizeId = getMipSizeId(m_blockSize);

  // init reduced boundary size
  m_reducedBdrySize = (m_sizeId == 0) ? 2 : 4;

  // init reduced prediction size
  m_reducedPredSize = (m_sizeId < 2) ? 4 : 8;

  // init upsampling factors
  m_upsmpFactorHor = m_blockSize.width / m_reducedPredSize;
  m_upsmpFactorVer = m_blockSize.height / m_reducedPredSize;

  CHECKD((m_upsmpFactorHor < 1) || ((m_upsmpFactorHor & (m_upsmpFactorHor - 1)) != 0),
         "Need power of two horizontal upsampling factor.");
  CHECKD((m_upsmpFactorVer < 1) || ((m_upsmpFactorVer & (m_upsmpFactorVer - 1)) != 0),
         "Need power of two vertical upsampling factor.");
}

void PredictorMIP::boundaryDownsampling1D(int16_t* reducedDst, const int16_t* const fullSrc, const SizeType srcLen,
                                          const SizeType dstLen) {
  if (dstLen < srcLen) {
    // Create reduced boundary by downsampling
    const SizeType downsmpFactor = srcLen / dstLen;
    const int log2DownsmpFactor = getLog2(downsmpFactor);
    const int roundingOffset = (1 << (log2DownsmpFactor - 1));

    SizeType srcIdx = 0;
    for (SizeType dstIdx = 0; dstIdx < dstLen; dstIdx++) {
      int sum = 0;
      for (int k = 0; k < downsmpFactor; k++) {
        sum += fullSrc[srcIdx++];
      }
      reducedDst[dstIdx] = (sum + roundingOffset) >> log2DownsmpFactor;
    }
  } else {
    // Copy boundary if no downsampling is needed
    for (SizeType i = 0; i < dstLen; ++i) {
      reducedDst[i] = fullSrc[i];
    }
  }
}

void PredictorMIP::predictionUpsampling1D(int16_t* const dst, const int16_t* const src, const int16_t* const bndry,
                                          const SizeType srcSizeUpsmpDim, const SizeType srcSizeOrthDim,
                                          const SizeType srcStep, const SizeType srcStride, const SizeType dstStep,
                                          const SizeType dstStride, const SizeType bndryStep,
                                          const unsigned int upsmpFactor) {
  const int log2UpsmpFactor = getLog2(upsmpFactor);
  CHECKD(upsmpFactor <= 1, "Upsampling factor must be at least 2.");
  const int roundingOffset = 1 << (log2UpsmpFactor - 1);

  SizeType idxOrthDim = 0;
  const int16_t* srcLine = src;
  int16_t* dstLine = dst;
  const int16_t* bndryLine = bndry + bndryStep - 1;
  while (idxOrthDim < srcSizeOrthDim) {
    SizeType idxUpsmpDim = 0;
    const int16_t* before = bndryLine;
    const int16_t* behind = srcLine;
    int16_t* currDst = dstLine;
    while (idxUpsmpDim < srcSizeUpsmpDim) {
      SizeType pos = 1;
      int scaledBefore = (*before) << log2UpsmpFactor;
      int scaledBehind = 0;
      while (pos <= upsmpFactor) {
        scaledBefore -= *before;
        scaledBehind += *behind;
        *currDst = (scaledBefore + scaledBehind + roundingOffset) >> log2UpsmpFactor;

        pos++;
        currDst += dstStep;
      }

      idxUpsmpDim++;
      before = behind;
      behind += srcStep;
    }

    idxOrthDim++;
    srcLine += srcStride;
    dstLine += dstStride;
    bndryLine += bndryStep;
  }
}

void PredictorMIP::predictionUpsampling(int16_t* const dst, const int16_t* const src) const {
  const int16_t* verSrc = src;
  SizeType verSrcStep = m_blockSize.width;

  if (m_upsmpFactorHor > 1) {
    int16_t* const horDst = dst + (m_upsmpFactorVer - 1) * m_blockSize.width;
    verSrc = horDst;
    verSrcStep *= m_upsmpFactorVer;

    upsampling1DHor(horDst, src, &m_refSamplesLeft[0], verSrcStep, m_upsmpFactorVer, m_reducedPredSize,
                    m_upsmpFactorHor);
  }

  if (m_upsmpFactorVer > 1) {
    upsampling1DVer(dst, verSrc, &m_refSamplesTop[0], m_blockSize.width, verSrcStep, m_reducedPredSize,
                    m_upsmpFactorVer);
  }
}

const uint8_t* PredictorMIP::getMatrixData(const int modeIdx) const {
  switch (m_sizeId) {
    case 0:
      return &mipMatrix4x4[modeIdx][0][0];

    case 1:
      return &mipMatrix8x8[modeIdx][0][0];

    case 2:
      return &mipMatrix16x16[modeIdx][0][0];

    default:
      THROW("Invalid mipSizeId");
  }
}

static void computeReducedPredCore(int16_t* const result, const int16_t* const input, const uint8_t* matrix,
                                   const bool transpose, const int bitDepth, const int inputSize, const int inputOffset,
                                   const int reducedPredSize, const int sizeId) {
  // use local buffer for transposed result
  int16_t resBufTransposed[MIP_MAX_REDUCED_OUTPUT_SAMPLES];
  int16_t* const resPtr = (transpose) ? resBufTransposed : result;

  int sum = 0;
  for (int i = 0; i < inputSize; i++) {
    sum += input[i];
  }
  const int offset = (1 << (MIP_SHIFT_MATRIX - 1)) - MIP_OFFSET_MATRIX * sum;
  CHECK(inputSize != 4 && inputSize != 8, "Error, inputSize must be 4 or 8");
  CHECK(reducedPredSize != 4 && reducedPredSize != 8, "Error, reducedPredSize must be 4 or 8");

  const uint8_t* weight = matrix;

  const bool redSize = (sizeId == 2);
  int posRes = 0;
  for (int y = 0; y < reducedPredSize; y++) {
    for (int x = 0; x < reducedPredSize; x++) {
      int tmp0 = redSize ? 0 : (input[0] * weight[0]);
      int tmp1 = input[1] * weight[1];
      int tmp2 = input[2] * weight[2];
      int tmp3 = input[3] * weight[3];
      for (int i = 4; i < inputSize; i += 4) {
        tmp0 += input[i] * weight[i];
        tmp1 += input[i + 1] * weight[i + 1];
        tmp2 += input[i + 2] * weight[i + 2];
        tmp3 += input[i + 3] * weight[i + 3];
      }
      resPtr[posRes++] =
          ClipBD<int>(((tmp0 + tmp1 + tmp2 + tmp3 + offset) >> MIP_SHIFT_MATRIX) + inputOffset, bitDepth);

      weight += inputSize;
    }
  }

  if (transpose) {
    for (int y = 0; y < reducedPredSize; y++) {
      for (int x = 0; x < reducedPredSize; x++) {
        result[y * reducedPredSize + x] = resPtr[x * reducedPredSize + y];
      }
    }
  }
}

static void upsampling1DVerCore(int16_t* const dst, const int16_t* const src, const int16_t* const bndry,
                                const int outWidth, const int srcStep, int inHeight, unsigned upsmpFactor) {
  const int roundingOffset = upsmpFactor >> 1;
  unsigned log2UpsmpFactor = getLog2(upsmpFactor);

  Pel* dstLine = dst;
  const Pel* srcLine = src;
  const Pel* bndryLine = bndry;

  for (int idxUpsmpDim = 0; idxUpsmpDim < inHeight; idxUpsmpDim++) {
    for (int idxOrthDim = 0; idxOrthDim < outWidth; idxOrthDim++) {
      const Pel* before = bndryLine + idxOrthDim;
      const Pel* behind = srcLine + idxOrthDim;
      Pel* currDst = dstLine + idxOrthDim;
      int scaledVal = (*before) << log2UpsmpFactor;
      for (int pos = 0; pos < upsmpFactor; pos++) {
        scaledVal -= *before;
        scaledVal += *behind;
        *currDst = (scaledVal + roundingOffset) >> log2UpsmpFactor;
        currDst += outWidth;
      }
    }
    dstLine += outWidth * upsmpFactor;
    bndryLine = srcLine;
    srcLine += srcStep;
  }
}

static void upsampling1DHorCore(int16_t* const dst, const int16_t* const src, const int16_t* const bndry,
                                const int dstStride, const int bndryStep, int predPredSize, unsigned upsmpFactor) {
  const int roundingOffset = upsmpFactor >> 1;
  unsigned log2UpsmpFactor = getLog2(upsmpFactor);

  Pel* dstLine = dst;
  const Pel* srcLine = src;
  const Pel* bndryLine = bndry + bndryStep - 1;

  for (SizeType idxOrthDim = 0; idxOrthDim < predPredSize; idxOrthDim++) {
    const Pel* before = bndryLine;
    const Pel* behind = srcLine;
    Pel* currDst = dstLine;
    for (SizeType idxUpsmpDim = 0; idxUpsmpDim < predPredSize; idxUpsmpDim++) {
      int scaledVal = (*before) << log2UpsmpFactor;
      for (SizeType pos = 0; pos < upsmpFactor; pos++) {
        scaledVal -= *before;
        scaledVal += *behind;
        *currDst = (scaledVal + roundingOffset) >> log2UpsmpFactor;
        currDst++;
      }
      before = behind;
      behind++;
    }

    srcLine += predPredSize;
    dstLine += dstStride;
    bndryLine += bndryStep;
  }
}

void PredictorMIP::initPredictorMIP(int betdepth) {
  computeReducedPred = computeReducedPredCore;
  upsampling1DVer = upsampling1DVerCore;
  upsampling1DHor = upsampling1DHorCore;
}
}  // namespace Mip

MatrixIntraPrediction::MatrixIntraPrediction() {}

void MatrixIntraPrediction::prepareInputForPred(const CPelBuf& src, const Area& puArea, const int bitDepth,
                                                const ComponentID compId) {
  m_component = compId;
  m_predictorMip.deriveBoundaryData(src, puArea, bitDepth);
}

void MatrixIntraPrediction::predBlock(const Size& puSize, const int intraMode, PelBuf& dst, const bool transpose,
                                      const int bitDepth, const ComponentID compId) {
  CHECK(m_component != compId, "Boundary has not been prepared for this component.");
  int16_t* const resultMip = m_mipResult;

  m_predictorMip.getPrediction(resultMip, intraMode, transpose, bitDepth);

#if ADAPTIVE_BIT_DEPTH
  if (bitDepth <= 8) {
    Pel8bit* dstPtr = reinterpret_cast<Pel8bit*>(dst.buf);
    ptrdiff_t dstStride = dst.stride;
    for (int y = 0; y < puSize.height; y++) {
      int16_t* const resultLine = &resultMip[y * puSize.width];
      Pel8bit* dstLine = dstPtr + y * dstStride;

      for (int x = 0; x < puSize.width; x += 4) {
        dstLine[x + 0] = Pel8bit(resultLine[x + 0]);
        dstLine[x + 1] = Pel8bit(resultLine[x + 1]);
        dstLine[x + 2] = Pel8bit(resultLine[x + 2]);
        dstLine[x + 3] = Pel8bit(resultLine[x + 3]);
      }
    }
  } else {
    for (int y = 0; y < puSize.height; y++) {
      int16_t* const resultLine = &resultMip[y * puSize.width];
      Pel* dstLine = dst.bufAt(0, y);

      for (int x = 0; x < puSize.width; x += 4) {
        dstLine[x + 0] = Pel(resultLine[x + 0]);
        dstLine[x + 1] = Pel(resultLine[x + 1]);
        dstLine[x + 2] = Pel(resultLine[x + 2]);
        dstLine[x + 3] = Pel(resultLine[x + 3]);
      }
    }
  }
#else
  for (int y = 0; y < puSize.height; y++) {
    int16_t* const resultLine = &resultMip[y * puSize.width];
    Pel* dstLine = dst.bufAt(0, y);

    for (int x = 0; x < puSize.width; x += 4) {
      dstLine[x + 0] = Pel(resultLine[x + 0]);
      dstLine[x + 1] = Pel(resultLine[x + 1]);
      dstLine[x + 2] = Pel(resultLine[x + 2]);
      dstLine[x + 3] = Pel(resultLine[x + 3]);
    }
  }
#endif
}
