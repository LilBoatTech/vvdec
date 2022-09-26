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

/** \file     SampleAdaptiveOffset.cpp
    \brief    sample adaptive offset class
*/

#include "SampleAdaptiveOffset.h"

#include "UnitTools.h"
#include "UnitPartitioner.h"
#include "CodingStructure.h"
#include "CommonLib/dtrace_codingstruct.h"
#include "CommonLib/dtrace_buffer.h"
#include "CommonLib/TimeProfiler.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//! \ingroup CommonLib
//! \{

template <typename TSrcDst>
static void processBO_c(Pel* _resLine, ptrdiff_t resStride, int width, int height, const int* offset, int /*startIdx*/,
                        const ClpRng& clpRng, int channelBitDepth) {
  TSrcDst* resLine = reinterpret_cast<TSrcDst*>(_resLine);
  const int shiftBits = channelBitDepth - NUM_SAO_BO_CLASSES_LOG2;
  int x, y;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      resLine[x] = ClipPel<int>(resLine[x] + offset[resLine[x] >> shiftBits], clpRng);
    }
    resLine += resStride;
  }
}

template <typename TSrcDst>
static void processEO0_c(Pel* _resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                         const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, const Pel* _leftLine) {
  TSrcDst* resLine = reinterpret_cast<TSrcDst*>(_resLine);
  const TSrcDst* leftLine = (const TSrcDst*)_leftLine;
  int startX = isLeftAvail ? 0 : 1;
  int endX = isRightAvail ? width : (width - 1);
  int signLeft;
  int x, y;
  offset += 2;
  for (y = 0; y < height; y++) {
    if (isLeftAvail)
      signLeft = (int8_t)sgn(resLine[startX] - leftLine[y]);
    else
      signLeft = (int8_t)sgn(resLine[startX] - resLine[startX - 1]);
    for (x = startX; x < endX; x++) {
      int signRight = (int8_t)sgn(resLine[x] - resLine[x + 1]);
      int edgeType = signRight + signLeft;
      signLeft = -signRight;
      resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
    }
    resLine += resStride;
  }
}

template <typename TSrcDst>
static void processEO90_c(Pel* _resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                          const ClpRng& clpRng, bool isAboveAvail, bool isBottomAvail, const Pel* _topLine) {
  TSrcDst* resLine = reinterpret_cast<TSrcDst*>(_resLine);
  const TSrcDst* topLine = (const TSrcDst*)_topLine;
  int8_t signUpLine[128];
  int startY = isAboveAvail ? 0 : 1;
  int endY = isBottomAvail ? height : (height - 1);
  int x, y;
  offset += 2;
  if (!isAboveAvail) {
    resLine += resStride;
  }
  const TSrcDst* srcLineAbove = resLine - resStride;
  if (isAboveAvail) srcLineAbove = topLine;
  for (x = 0; x < width; x++) {
    signUpLine[x] = (int8_t)sgn(resLine[x] - srcLineAbove[x]);
  }
  for (y = startY; y < endY; y++) {
    const TSrcDst* srcLineBelow = resLine + resStride;
    for (x = 0; x < width; x++) {
      int8_t signDown = (int8_t)sgn(resLine[x] - srcLineBelow[x]);
      int edgeType = signDown + signUpLine[x];
      signUpLine[x] = -signDown;
      resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
    }
    resLine += resStride;
  }
}

template <typename TSrcDst>
static void processEO135_c(Pel* _resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                           const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, bool isAboveLeftAvail,
                           bool isAboveAvail, bool isBelowAvail, bool isBelowRightAvail, const Pel* _topLine,
                           const Pel* _leftLine) {
  TSrcDst* resLine = reinterpret_cast<TSrcDst*>(_resLine);
  const TSrcDst* topLine = (const TSrcDst*)_topLine;
  const TSrcDst* leftLine = (const TSrcDst*)_leftLine;
  offset += 2;
  int8_t *signUpLine, *signDownLine, *signTmpLine;
  int8_t aSignUpLine[128 + 1], aSignDownLine[128 + 1];
  signUpLine = &aSignUpLine[0];
  signDownLine = &aSignDownLine[0];

  int startX = isLeftAvail ? 0 : 1;
  int endX = isRightAvail ? width : (width - 1);
  int x, y;
  // prepare 2nd line's upper sign
  const TSrcDst* srcLineBelow = resLine + resStride;
  x = startX;
  if (isLeftAvail) {
    signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - leftLine[0]);
    x++;
  }
  for (; x < endX + 1; x++) {
    signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - resLine[x - 1]);
  }

  // 1st line
  const TSrcDst* srcLineAbove = resLine - resStride;
  int firstLineStartX = isAboveLeftAvail ? 0 : 1;
  int firstLineEndX = isAboveAvail ? endX : 1;
  if (isAboveAvail) srcLineAbove = topLine;
  for (x = firstLineStartX; x < firstLineEndX; x++) {
    int edgeType = sgn(resLine[x] - srcLineAbove[x - 1]) - signUpLine[x + 1];
    resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
  }
  resLine += resStride;
  // middle lines
  for (y = 1; y < height - 1; y++) {
    srcLineBelow = resLine + resStride;
    for (x = startX; x < endX; x++) {
      int signDown = (int8_t)sgn(resLine[x] - srcLineBelow[x + 1]);
      int edgeType = signDown + signUpLine[x];
      resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);

      signDownLine[x + 1] = -signDown;
    }
    if (isLeftAvail)
      signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - leftLine[y]);
    else
      signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - resLine[startX - 1]);
    signTmpLine = signUpLine;
    signUpLine = signDownLine;
    signDownLine = signTmpLine;
    resLine += resStride;
  }

  // last line
  srcLineBelow = resLine + resStride;
  int lastLineStartX = isBelowAvail ? startX : (width - 1);
  int lastLineEndX = isBelowRightAvail ? width : (width - 1);
  for (x = lastLineStartX; x < lastLineEndX; x++) {
    int edgeType = sgn(resLine[x] - srcLineBelow[x + 1]) + signUpLine[x];
    resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
  }
}

template <typename TSrcDst>
static void processEO45_c(Pel* _resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                          const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, bool isAboveRightAvail,
                          bool isAboveAvail, bool isBelowAvail, bool isBelowLeftAvail, const Pel* _topLine,
                          const Pel* _leftLine) {
  TSrcDst* resLine = reinterpret_cast<TSrcDst*>(_resLine);
  const TSrcDst* topLine = (const TSrcDst*)_topLine;
  const TSrcDst* leftLine = (const TSrcDst*)_leftLine;
  offset += 2;
  int8_t aSignUpLine[128 + 2];
  int8_t* signUpLine = &aSignUpLine[1];
  int x, y;
  int startX = isLeftAvail ? 0 : 1;
  int endX = isRightAvail ? width : (width - 1);

  // prepare 2nd line upper sign
  const TSrcDst* srcLineBelow = resLine + resStride;
  x = startX - 1;
  if (isLeftAvail) {
    signUpLine[x] = (int8_t)sgn(leftLine[1] - resLine[x + 1]);
    x++;
  }
  for (; x < endX; x++) {
    signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - resLine[x + 1]);
  }
  // first line
  const TSrcDst* srcLineAbove = resLine - resStride;
  int firstLineStartX = isAboveAvail ? startX : (width - 1);
  int firstLineEndX = isAboveRightAvail ? width : (width - 1);
  if (isAboveAvail) srcLineAbove = topLine;
  for (x = firstLineStartX; x < firstLineEndX; x++) {
    int edgeType = sgn(resLine[x] - srcLineAbove[x + 1]) - signUpLine[x - 1];
    resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
  }
  resLine += resStride;

  // middle lines
  for (y = 1; y < height - 1; y++) {
    srcLineBelow = resLine + resStride;
    x = startX;
    if (isLeftAvail) {
      int signDown = (int8_t)sgn(resLine[x] - leftLine[y + 1]);
      int edgeType = signDown + signUpLine[x];
      resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      signUpLine[x - 1] = -signDown;
      x++;
    }
    for (; x < endX; x++) {
      int signDown = (int8_t)sgn(resLine[x] - srcLineBelow[x - 1]);
      int edgeType = signDown + signUpLine[x];
      resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      signUpLine[x - 1] = -signDown;
    }
    signUpLine[endX - 1] = (int8_t)sgn(srcLineBelow[endX - 1] - resLine[endX]);
    resLine += resStride;
  }

  // last line
  srcLineBelow = resLine + resStride;
  int lastLineStartX = isBelowLeftAvail ? 0 : 1;
  int lastLineEndX = isBelowAvail ? endX : 1;
  // below left is in reconstruction buffer, not line buffer
  for (x = lastLineStartX; x < lastLineEndX; x++) {
    int edgeType = sgn(resLine[x] - srcLineBelow[x - 1]) + signUpLine[x];
    resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
  }
}

void SampleAdaptiveOffset::offsetBlock_core(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset,
                                            int startIdx, Pel* resBlk, ptrdiff_t resStride, int width, int height,
                                            bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail,
                                            bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail,
                                            bool isBelowRightAvail, const Pel* topline, const Pel* leftline) {
  switch (typeIdx) {
    case SAO_TYPE_EO_0: {
      processEO0(resBlk, resStride, width, height, offset, clpRng, isLeftAvail, isRightAvail, leftline);
      break;
    }
    case SAO_TYPE_EO_90: {
      processEO90(resBlk, resStride, width, height, offset, clpRng, isAboveAvail, isBelowAvail, topline);
      break;
    }
    case SAO_TYPE_EO_135: {
      processEO135(resBlk, resStride, width, height, offset, clpRng, isLeftAvail, isRightAvail, isAboveLeftAvail,
                   isAboveAvail, isBelowAvail, isBelowRightAvail, topline, leftline);
      break;
    }
    case SAO_TYPE_EO_45: {
      processEO45(resBlk, resStride, width, height, offset, clpRng, isLeftAvail, isRightAvail, isAboveRightAvail,
                  isAboveAvail, isBelowAvail, isBelowLeftAvail, topline, leftline);
      break;
    }
    case SAO_TYPE_BO: {
      processBO(resBlk, resStride, width, height, offset, startIdx, clpRng, channelBitDepth);
      break;
    }
  }
}

template <typename T>
void SampleAdaptiveOffset::offsetBlock_core(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset,
                                            int startIdx, T* resBlk, ptrdiff_t resStride, int width, int height,
                                            bool isLeftAvail, bool isRightAvail, bool isAboveAvail, bool isBelowAvail,
                                            bool isAboveLeftAvail, bool isAboveRightAvail, bool isBelowLeftAvail,
                                            bool isBelowRightAvail, bool isCtuCrossedByVirtualBoundaries,
                                            int horVirBndryPos[], int verVirBndryPos[], int numHorVirBndry,
                                            int numVerVirBndry, const T* topline, const T* leftline) {
  int x, y, startX, startY, endX, endY, edgeType;
  int firstLineStartX, firstLineEndX, lastLineStartX, lastLineEndX;
  int8_t signLeft, signRight, signDown;
  T* resLine = reinterpret_cast<T*>(resBlk);

  switch (typeIdx) {
    case SAO_TYPE_EO_0: {
      offset += 2;
      startX = isLeftAvail ? 0 : 1;
      endX = isRightAvail ? width : (width - 1);
      for (y = 0; y < height; y++) {
        if (isLeftAvail)
          signLeft = (int8_t)sgn(resLine[startX] - leftline[y]);
        else
          signLeft = (int8_t)sgn(resLine[startX] - resLine[startX - 1]);
        for (x = startX; x < endX; x++) {
          signRight = (int8_t)sgn(resLine[x] - resLine[x + 1]);
          if (isCtuCrossedByVirtualBoundaries &&
              isProcessDisabled(x, y, numVerVirBndry, 0, verVirBndryPos, horVirBndryPos)) {
            signLeft = -signRight;
            continue;
          }
          edgeType = signRight + signLeft;
          signLeft = -signRight;

          resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
        }
        resLine += resStride;
      }
    } break;

    case SAO_TYPE_EO_90: {
      offset += 2;
      int8_t aSignUpLine[128];
      int8_t* signUpLine = &aSignUpLine[0];
      startY = isAboveAvail ? 0 : 1;
      endY = isBelowAvail ? height : height - 1;
      if (!isAboveAvail) {
        resLine += resStride;
      }

      const T* srcLineAbove = resLine - resStride;
      if (isAboveAvail) srcLineAbove = topline;
      for (x = 0; x < width; x++) {
        signUpLine[x] = (int8_t)sgn(resLine[x] - srcLineAbove[x]);
      }

      const T* srcLineBelow;
      for (y = startY; y < endY; y++) {
        srcLineBelow = resLine + resStride;

        for (x = 0; x < width; x++) {
          signDown = (int8_t)sgn(resLine[x] - srcLineBelow[x]);
          if (isCtuCrossedByVirtualBoundaries &&
              isProcessDisabled(x, y, 0, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
            signUpLine[x] = -signDown;
            continue;
          }
          edgeType = signDown + signUpLine[x];
          signUpLine[x] = -signDown;

          resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
        }
        resLine += resStride;
      }
    } break;
    case SAO_TYPE_EO_135: {
      offset += 2;
      int8_t *signUpLine, *signDownLine, *signTmpLine;
      int8_t aSignUpLine[128 + 1], aSignDownLine[128 + 1];
      signUpLine = &aSignUpLine[0];
      signDownLine = &aSignDownLine[0];

      startX = isLeftAvail ? 0 : 1;
      endX = isRightAvail ? width : (width - 1);

      // prepare 2nd line's upper sign
      const T* srcLineBelow = resLine + resStride;
      x = startX;
      if (isLeftAvail) {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - leftline[0]);
        x++;
      }
      for (; x < endX + 1; x++) {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - resLine[x - 1]);
      }

      // 1st line
      const T* srcLineAbove = resLine - resStride;
      firstLineStartX = isAboveLeftAvail ? 0 : 1;
      firstLineEndX = isAboveAvail ? endX : 1;
      if (isAboveAvail) srcLineAbove = topline;
      for (x = firstLineStartX; x < firstLineEndX; x++) {
        if (isCtuCrossedByVirtualBoundaries &&
            isProcessDisabled(x, 0, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
          continue;
        }
        edgeType = sgn(resLine[x] - srcLineAbove[x - 1]) - signUpLine[x + 1];

        resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      }
      resLine += resStride;

      // middle lines
      for (y = 1; y < height - 1; y++) {
        srcLineBelow = resLine + resStride;

        for (x = startX; x < endX; x++) {
          signDown = (int8_t)sgn(resLine[x] - srcLineBelow[x + 1]);
          if (isCtuCrossedByVirtualBoundaries &&
              isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
            signDownLine[x + 1] = -signDown;
            continue;
          }
          edgeType = signDown + signUpLine[x];
          resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);

          signDownLine[x + 1] = -signDown;
        }
        if (isLeftAvail)
          signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - leftline[y]);
        else
          signDownLine[startX] = (int8_t)sgn(srcLineBelow[startX] - resLine[startX - 1]);

        signTmpLine = signUpLine;
        signUpLine = signDownLine;
        signDownLine = signTmpLine;

        resLine += resStride;
      }

      // last line
      srcLineBelow = resLine + resStride;
      lastLineStartX = isBelowAvail ? startX : (width - 1);
      lastLineEndX = isBelowRightAvail ? width : (width - 1);
      for (x = lastLineStartX; x < lastLineEndX; x++) {
        if (isCtuCrossedByVirtualBoundaries &&
            isProcessDisabled(x, height - 1, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
          continue;
        }
        edgeType = sgn(resLine[x] - srcLineBelow[x + 1]) + signUpLine[x];
        resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      }
    } break;
    case SAO_TYPE_EO_45: {
      offset += 2;
      int8_t aSignUpLine[128 + 2];
      int8_t* signUpLine = &aSignUpLine[1];

      startX = isLeftAvail ? 0 : 1;
      endX = isRightAvail ? width : (width - 1);

      // prepare 2nd line upper sign
      const T* srcLineBelow = resLine + resStride;
      x = startX - 1;
      if (isLeftAvail) {
        signUpLine[x] = (int8_t)sgn(leftline[1] - resLine[x + 1]);
        x++;
      }
      for (; x < endX; x++) {
        signUpLine[x] = (int8_t)sgn(srcLineBelow[x] - resLine[x + 1]);
      }

      // first line
      const T* srcLineAbove = resLine - resStride;
      firstLineStartX = isAboveAvail ? startX : (width - 1);
      firstLineEndX = isAboveRightAvail ? width : (width - 1);
      if (isAboveAvail) srcLineAbove = topline;
      for (x = firstLineStartX; x < firstLineEndX; x++) {
        if (isCtuCrossedByVirtualBoundaries &&
            isProcessDisabled(x, 0, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
          continue;
        }
        edgeType = sgn(resLine[x] - srcLineAbove[x + 1]) - signUpLine[x - 1];
        resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      }
      resLine += resStride;

      // middle lines
      for (y = 1; y < height - 1; y++) {
        srcLineBelow = resLine + resStride;
        x = startX;
        if (isLeftAvail) {
          signDown = (int8_t)sgn(resLine[x] - leftline[y + 1]);
          if (isCtuCrossedByVirtualBoundaries &&
              isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
            signUpLine[x - 1] = -signDown;
          } else {
            edgeType = signDown + signUpLine[x];
            resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
            signUpLine[x - 1] = -signDown;
            x++;
          }
        }
        for (; x < endX; x++) {
          signDown = (int8_t)sgn(resLine[x] - srcLineBelow[x - 1]);
          if (isCtuCrossedByVirtualBoundaries &&
              isProcessDisabled(x, y, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
            signUpLine[x - 1] = -signDown;
            continue;
          }
          edgeType = signDown + signUpLine[x];
          resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
          signUpLine[x - 1] = -signDown;
        }
        signUpLine[endX - 1] = (int8_t)sgn(srcLineBelow[endX - 1] - resLine[endX]);
        resLine += resStride;
      }

      // last line
      srcLineBelow = resLine + resStride;
      lastLineStartX = isBelowLeftAvail ? 0 : 1;
      lastLineEndX = isBelowAvail ? endX : 1;
      // below left is in reconstruction buffer, not line buffer
      for (x = lastLineStartX; x < lastLineEndX; x++) {
        if (isCtuCrossedByVirtualBoundaries &&
            isProcessDisabled(x, height - 1, numVerVirBndry, numHorVirBndry, verVirBndryPos, horVirBndryPos)) {
          continue;
        }
        edgeType = sgn(resLine[x] - srcLineBelow[x - 1]) + signUpLine[x];
        resLine[x] = ClipPel<int>(resLine[x] + offset[edgeType], clpRng);
      }
    } break;
    case SAO_TYPE_BO: {
      const int shiftBits = channelBitDepth - NUM_SAO_BO_CLASSES_LOG2;

      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          resLine[x] = ClipPel<int>(resLine[x] + offset[resLine[x] >> shiftBits], clpRng);
        }
        resLine += resStride;
      }
    } break;
    default: {
      // THROW("Not a supported SAO types\n");
    }
  }
}

SampleAdaptiveOffset::~SampleAdaptiveOffset() { destroy(); }

void SampleAdaptiveOffset::initSampleAdaptiveOffset(int bitDepth) {
#define INIT_SAO_FUNC(type)            \
  processBO = processBO_c<type>;       \
  processEO0 = processEO0_c<type>;     \
  processEO90 = processEO90_c<type>;   \
  processEO135 = processEO135_c<type>; \
  processEO45 = processEO45_c<type>;

#if ADAPTIVE_BIT_DEPTH
  if (bitDepth <= 8) {
    INIT_SAO_FUNC(Pel8bit);
  } else {
    INIT_SAO_FUNC(Pel);
  }
#else
  INIT_SAO_FUNC(Pel);
#endif

#undef INIT_SAO_FUNC

  //  processBO = processBO_c<Pel>;  // c code
  //  processEO0 = processEO0_c<Pel>;
  //  processEO90 = processEO90_c<Pel>;
  //  processEO135 = processEO135_c<Pel>;
  //  processEO45 = processEO45_c<Pel>;
}

void SampleAdaptiveOffset::create(int picWidth, int picHeight, ChromaFormat format, uint32_t maxCUWidth,
                                  uint32_t maxCUHeight, uint32_t maxCUDepth) {
  if (picWidth != m_picW || picHeight != m_picH || format != m_chromaFormat)
    destroy();
  else
    return;
  m_picW = picWidth;
  m_picH = picHeight;
  m_chromaFormat = format;

  m_numberOfComponents = getNumberValidComponents(format);

  int numCTUCol = (picWidth + maxCUWidth - 1) / maxCUWidth;
  int numCTURow = (picHeight + maxCUHeight - 1) / maxCUHeight;
  for (int compIdx = 0; compIdx < m_numberOfComponents; compIdx++) {
    m_lineTopBuf[compIdx] = xMalloc(Pel*, numCTURow);
    m_lineLeftBuf[compIdx] = xMalloc(Pel*, numCTUCol);
    int rowLength = (picWidth >> getComponentScaleX(ComponentID(compIdx), format)) + 16;
    m_lineTopBuf[compIdx][0] = xMalloc(Pel, rowLength * numCTURow);
    for (int r = 1; r < numCTURow; r++) m_lineTopBuf[compIdx][r] = m_lineTopBuf[compIdx][r - 1] + rowLength;
    int colLength = (picHeight >> getComponentScaleY(ComponentID(compIdx), format)) + 16;
    m_lineLeftBuf[compIdx][0] = xMalloc(Pel, colLength * numCTUCol);
    for (int c = 1; c < numCTUCol; c++) m_lineLeftBuf[compIdx][c] = m_lineLeftBuf[compIdx][c - 1] + colLength;
  }
}

void SampleAdaptiveOffset::destroy() {
  for (int compIdx = 0; compIdx < m_numberOfComponents; compIdx++) {
    if (m_lineTopBuf[compIdx] != nullptr) xFree(m_lineTopBuf[compIdx][0]);
    if (m_lineLeftBuf[compIdx] != nullptr) xFree(m_lineLeftBuf[compIdx][0]);
    xFree(m_lineTopBuf[compIdx]);
    xFree(m_lineLeftBuf[compIdx]);
    m_lineTopBuf[compIdx] = nullptr;
    m_lineLeftBuf[compIdx] = nullptr;
  }
}

void SampleAdaptiveOffset::SAOPrepareCTULine(CodingStructure& cs, const UnitArea& lineArea) {
  PROFILER_SCOPE_AND_STAGE(1, g_timeProfiler, P_SAO);

  const PreCalcValues& pcv = *cs.pcv;
  const int y = lineArea.lumaPos().y;
  for (int x = 0; x < pcv.lumaWidth; x += pcv.maxCUWidth) {
    const int ctuRsAddr = getCtuAddr(Position(x, y), *cs.pcv);

    SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES] = {nullptr, nullptr};
    getMergeList(cs, ctuRsAddr, mergeList);
    reconstructBlkSAOParam(cs.m_saoBlkParam[ctuRsAddr], mergeList);
  }
}

void SampleAdaptiveOffset::SAOProcessCTU(CodingStructure& cs, const UnitArea& ctuArea) {
  PROFILER_SCOPE_AND_STAGE(1, g_timeProfiler, P_SAO);
  const int ctuRsAddr = getCtuAddr(ctuArea.lumaPos(), *cs.pcv);
  // SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES] = { nullptr, nullptr };
  // getMergeList( cs, ctuRsAddr, mergeList );
  // reconstructBlkSAOParam( cs.m_saoBlkParam[ctuRsAddr], mergeList );
  offsetCTU(ctuArea, *cs.getReconBufPtr(), ctuRsAddr, cs);
}

int SampleAdaptiveOffset::getMergeList(CodingStructure& cs, int ctuRsAddr,
                                       SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES]) {
  const PreCalcValues& pcv = *cs.pcv;

  const int ctuX = ctuRsAddr % pcv.widthInCtus;
  const int ctuY = ctuRsAddr / pcv.widthInCtus;
  const CodingUnit& cu = *cs.getCtuCuPtrData(ctuRsAddr).cuPtr[CH_L][0];
  int mergedCTUPos;
  int numValidMergeCandidates = 0;

  for (int mergeType = 0; mergeType < NUM_SAO_MERGE_TYPES; mergeType++) {
    SAOBlkParam* mergeCandidate = NULL;

    switch (mergeType) {
      case SAO_MERGE_ABOVE: {
        if (ctuY > 0) {
          mergedCTUPos = ctuRsAddr - pcv.widthInCtus;
          if (cu.above) {
            mergeCandidate = &(cs.m_saoBlkParam[mergedCTUPos]);
          }
        }
      } break;
      case SAO_MERGE_LEFT: {
        if (ctuX > 0) {
          mergedCTUPos = ctuRsAddr - 1;
          if (cu.left) {
            mergeCandidate = &(cs.m_saoBlkParam[mergedCTUPos]);
          }
        }
      } break;
      default: {
        THROW("not a supported merge type");
      }
    }

    mergeList[mergeType] = mergeCandidate;
    if (mergeCandidate != NULL) {
      numValidMergeCandidates++;
    }
  }

  return numValidMergeCandidates;
}

void SampleAdaptiveOffset::reconstructBlkSAOParam(SAOBlkParam& recParam,
                                                  SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES]) const {
  const int numberOfComponents = m_numberOfComponents;
  for (int compIdx = 0; compIdx < numberOfComponents; compIdx++) {
    const ComponentID component = ComponentID(compIdx);
    SAOOffset& offsetParam = recParam[component];

    if (offsetParam.modeIdc == SAO_MODE_OFF) {
      continue;
    }

    switch (offsetParam.modeIdc) {
      case SAO_MODE_NEW: {
        // invertQuantOffsets(component, offsetParam.typeIdc, offsetParam.typeAuxInfo, offsetParam.offset,
        // offsetParam.offset);
      } break;
      case SAO_MODE_MERGE: {
        const SAOBlkParam* mergeTarget = mergeList[offsetParam.typeIdc];
        CHECK(mergeTarget == NULL, "Merge target does not exist");

        offsetParam = (*mergeTarget)[component];
      } break;
      default: {
        THROW("Not a supported mode");
      }
    }
  }
}

template <typename T>
void SampleAdaptiveOffset::saveCTULine(CodingStructure& cs, const UnitArea& area, int compIdx, bool bSaveBottom,
                                       bool bSaveRight, bool bSaveRightFast) {
#if ADAPTIVE_BIT_DEPTH
  int bytePerPixel = cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) <= 8 ? 1 : 2;
#endif
  const ComponentID compID = ComponentID(compIdx);
  const CompArea& compArea = area.block(compID);
  PelBuf srcPlane = cs.getRecoBuf(compID);
  int width = compArea.width;
  int height = compArea.height;
  ptrdiff_t stride = srcPlane.stride;
  if (bSaveBottom) {
    int ctuRowIdx = (area.ly() >> cs.pcv->maxCUHeightLog2) + 1;
    T* dst = (reinterpret_cast<T*>(m_lineTopBuf[compIdx][ctuRowIdx])) + compArea.pos().x;
#if ADAPTIVE_BIT_DEPTH
    const T* src = (reinterpret_cast<T*>(srcPlane.bufAt(compArea.pos(), bytePerPixel))) + stride * (height - 1);
#else
    const T* src = (reinterpret_cast<T*>(srcPlane.bufAt(compArea.pos()))) + stride * (height - 1);
#endif
    memcpy(dst, src, width * sizeof(T));
  }
  if (bSaveRight) {
    int ctuColIdx = (area.lx() >> cs.pcv->maxCUWidthLog2) + 1;
    T* dst = (reinterpret_cast<T*>(m_lineLeftBuf[compIdx][ctuColIdx])) + compArea.pos().y;
#if ADAPTIVE_BIT_DEPTH
    const T* src = (reinterpret_cast<T*>(srcPlane.bufAt(compArea.pos(), bytePerPixel))) + (width - 1);
#else
    const T* src = (reinterpret_cast<T*>(srcPlane.bufAt(compArea.pos()))) + (width - 1);
#endif
    if (!bSaveRightFast) {
      for (int i = 0; i < height; i++) {
        dst[i] = src[0];
        src += stride;
      }
    } else {
      // only save 2 points
      dst[0] = src[0];
      dst[height - 1] = src[(height - 1) * stride];
    }
  }
}

void SampleAdaptiveOffset::offsetCTU(const UnitArea& area, PelUnitBuf& res, int ctuRSAddr, CodingStructure& cs) {
  const uint32_t numberOfComponents = getNumberValidComponents(area.chromaFormat);

#if ADAPTIVE_BIT_DEPTH
  int bytePerPixel = cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) <= 8 ? 1 : 2;
#endif

  bool isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail,
      isBelowRightAvail;
  SAOBlkParam& saoblkParam = cs.m_saoBlkParam[ctuRSAddr];
  // block boundary availability
  deriveLoopFilterBoundaryAvailibility(cs, area.Y(), isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail,
                                       isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail);

  int numHorVirBndry = 0;
  int numVerVirBndry = 0;
  int horVirBndryPos[] = {-1, -1, -1};
  int verVirBndryPos[] = {-1, -1, -1};
  int horVirBndryPosComp[] = {-1, -1, -1};
  int verVirBndryPosComp[] = {-1, -1, -1};
  bool isCtuCrossedByVirtualBoundaries = isCrossedByVirtualBoundaries(cs.picHeader, area.Y(), numHorVirBndry,
                                                                      numVerVirBndry, horVirBndryPos, verVirBndryPos);
  if (isCtuCrossedByVirtualBoundaries) {
    CHECK(numHorVirBndry >= (int)(sizeof(horVirBndryPos) / sizeof(horVirBndryPos[0])), "Too many virtual boundaries");
    CHECK(numHorVirBndry >= (int)(sizeof(verVirBndryPos) / sizeof(verVirBndryPos[0])), "Too many virtual boundaries");
  }

  int rowIdx = area.ly() >> cs.pcv->maxCUHeightLog2;
  int colIdx = area.lx() >> cs.pcv->maxCUWidthLog2;

  for (int compIdx = 0; compIdx < numberOfComponents; compIdx++) {
    const ComponentID compID = ComponentID(compIdx);
    const CompArea& compArea = area.block(compID);
    SAOOffset& ctbOffset = saoblkParam[compIdx];
    ptrdiff_t resStride = res.get(compID).stride;
#if ADAPTIVE_BIT_DEPTH
    Pel* resBlk = res.get(compID).bufAt(compArea, bytePerPixel);
#else
    Pel* resBlk = res.get(compID).bufAt(compArea);
#endif
    // save ctu line
    bool bSaveRightFast = false;
    if (isRightAvail) {
      const SAOOffset& next = cs.m_saoBlkParam[ctuRSAddr + 1][compIdx];
      if (next.modeIdc == SAO_MODE_OFF || next.typeIdc == SAO_TYPE_BO || next.typeIdc == SAO_TYPE_EO_90)
        bSaveRightFast = true;
    }
#if ADAPTIVE_BIT_DEPTH
    if (bytePerPixel == 1) {
      saveCTULine<Pel8bit>(cs, area, compIdx, isBelowAvail, isRightAvail, bSaveRightFast);
    } else {
      saveCTULine<Pel>(cs, area, compIdx, isBelowAvail, isRightAvail, bSaveRightFast);
    }
#else
    saveCTULine<Pel>(cs, area, compIdx, isBelowAvail, isRightAvail, bSaveRightFast);
#endif
    for (int i = 0; i < numHorVirBndry; i++) {
      horVirBndryPosComp[i] = (horVirBndryPos[i] >> ::getComponentScaleY(compID, area.chromaFormat)) - compArea.y;
    }
    for (int i = 0; i < numVerVirBndry; i++) {
      verVirBndryPosComp[i] = (verVirBndryPos[i] >> ::getComponentScaleX(compID, area.chromaFormat)) - compArea.x;
    }
    if (ctbOffset.modeIdc != SAO_MODE_OFF) {
      if (isCtuCrossedByVirtualBoundaries) {
#if ADAPTIVE_BIT_DEPTH
        if (bytePerPixel == 1) {
          offsetBlock_core<Pel8bit>(
              cs.sps->getBitDepth(toChannelType(compID)),
              cs.getCtuCuPtrData(cs.ctuRsAddr(area.Y().pos(), CH_L)).cuPtr[0][0]->slice->clpRng(compID),
              ctbOffset.typeIdc, ctbOffset.offset, ctbOffset.typeAuxInfo, reinterpret_cast<Pel8bit*>(resBlk), resStride,
              compArea.width, compArea.height, isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail,
              isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail, isCtuCrossedByVirtualBoundaries,
              horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry,
              (reinterpret_cast<Pel8bit*>(m_lineTopBuf[compIdx][rowIdx])) + compArea.pos().x,
              (reinterpret_cast<Pel8bit*>(m_lineLeftBuf[compIdx][colIdx])) + compArea.pos().y);
        } else {
          offsetBlock_core<Pel>(
              cs.sps->getBitDepth(toChannelType(compID)),
              cs.getCtuCuPtrData(cs.ctuRsAddr(area.Y().pos(), CH_L)).cuPtr[0][0]->slice->clpRng(compID),
              ctbOffset.typeIdc, ctbOffset.offset, ctbOffset.typeAuxInfo, resBlk, resStride, compArea.width,
              compArea.height, isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail,
              isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail, isCtuCrossedByVirtualBoundaries,
              horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry, numVerVirBndry,
              m_lineTopBuf[compIdx][rowIdx] + compArea.pos().x, m_lineLeftBuf[compIdx][colIdx] + compArea.pos().y);
        }
#else
        offsetBlock_core<Pel>(cs.sps->getBitDepth(toChannelType(compID)),
                              cs.getCtuCuPtrData(cs.ctuRsAddr(area.Y().pos(), CH_L)).cuPtr[0][0]->slice->clpRng(compID),
                              ctbOffset.typeIdc, ctbOffset.offset, ctbOffset.typeAuxInfo, resBlk, resStride,
                              compArea.width, compArea.height, isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail,
                              isAboveLeftAvail, isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail,
                              isCtuCrossedByVirtualBoundaries, horVirBndryPosComp, verVirBndryPosComp, numHorVirBndry,
                              numVerVirBndry, m_lineTopBuf[compIdx][rowIdx] + compArea.pos().x,
                              m_lineLeftBuf[compIdx][colIdx] + compArea.pos().y);
#endif
      } else {  // fast path
#if ADAPTIVE_BIT_DEPTH
        const Pel* topline =
            (const Pel*)((reinterpret_cast<uint8_t*>(m_lineTopBuf[compIdx][rowIdx])) + compArea.pos().x * bytePerPixel);
        const Pel* leftline = (const Pel*)((reinterpret_cast<uint8_t*>(m_lineLeftBuf[compIdx][colIdx])) +
                                           compArea.pos().y * bytePerPixel);
#else
        const Pel* topline = m_lineTopBuf[compIdx][rowIdx] + compArea.pos().x;
        const Pel* leftline = m_lineLeftBuf[compIdx][colIdx] + compArea.pos().y;
#endif
        offsetBlock_core(cs.sps->getBitDepth(toChannelType(compID)),
                         cs.getCtuCuPtrData(cs.ctuRsAddr(area.Y().pos(), CH_L)).cuPtr[0][0]->slice->clpRng(compID),
                         ctbOffset.typeIdc, ctbOffset.offset, ctbOffset.typeAuxInfo, resBlk, resStride, compArea.width,
                         compArea.height, isLeftAvail, isRightAvail, isAboveAvail, isBelowAvail, isAboveLeftAvail,
                         isAboveRightAvail, isBelowLeftAvail, isBelowRightAvail, topline, leftline);
      }
    }
  }
}

void SampleAdaptiveOffset::deriveLoopFilterBoundaryAvailibility(CodingStructure& cs, const Position& pos,
                                                                bool& isLeftAvail, bool& isRightAvail,
                                                                bool& isAboveAvail, bool& isBelowAvail,
                                                                bool& isAboveLeftAvail, bool& isAboveRightAvail,
                                                                bool& isBelowLeftAvail, bool& isBelowRightAvail) const {
  const int ctusz = cs.pcv->maxCUWidth;
  const int ctuX = pos.x / ctusz;
  const int ctuY = pos.y / ctusz;
  const int width = cs.pcv->widthInCtus;
  const int height = cs.pcv->heightInCtus;

  const CodingUnit* cuCurr = cs.getCtuCuPtrData(ctuX, ctuY).cuPtr[0][0];
  const CodingUnit* cuLeft = ctuX > 0 ? cs.getCtuCuPtrData(ctuX - 1, ctuY).cuPtr[0][0] : nullptr;
  const CodingUnit* cuRight = ctuX + 1 < width ? cs.getCtuCuPtrData(ctuX + 1, ctuY).cuPtr[0][0] : nullptr;
  const CodingUnit* cuAbove = ctuY > 0 ? cs.getCtuCuPtrData(ctuX, ctuY - 1).cuPtr[0][0] : nullptr;
  const CodingUnit* cuBelow = ctuY + 1 < height ? cs.getCtuCuPtrData(ctuX, ctuY + 1).cuPtr[0][0] : nullptr;
  const CodingUnit* cuAboveLeft = cuLeft && cuAbove ? cs.getCtuCuPtrData(ctuX - 1, ctuY - 1).cuPtr[0][0] : nullptr;
  const CodingUnit* cuAboveRight = cuRight && cuAbove ? cs.getCtuCuPtrData(ctuX + 1, ctuY - 1).cuPtr[0][0] : nullptr;
  const CodingUnit* cuBelowLeft = cuLeft && cuBelow ? cs.getCtuCuPtrData(ctuX - 1, ctuY + 1).cuPtr[0][0] : nullptr;
  const CodingUnit* cuBelowRight = cuRight && cuBelow ? cs.getCtuCuPtrData(ctuX + 1, ctuY + 1).cuPtr[0][0] : nullptr;

  // check cross slice flags
  const bool isLoopFilterAcrossSlicePPS = cs.pps->getLoopFilterAcrossSlicesEnabledFlag();
  if (!isLoopFilterAcrossSlicePPS) {
    isLeftAvail = (cuLeft == NULL) ? false : CU::isSameSlice(*cuCurr, *cuLeft);
    isAboveAvail = (cuAbove == NULL) ? false : CU::isSameSlice(*cuCurr, *cuAbove);
    isRightAvail = (cuRight == NULL) ? false : CU::isSameSlice(*cuCurr, *cuRight);
    isBelowAvail = (cuBelow == NULL) ? false : CU::isSameSlice(*cuCurr, *cuBelow);
    isAboveLeftAvail = (cuAboveLeft == NULL) ? false : CU::isSameSlice(*cuCurr, *cuAboveLeft);
    isAboveRightAvail = (cuAboveRight == NULL) ? false : CU::isSameSlice(*cuCurr, *cuAboveRight);
    isBelowLeftAvail = (cuBelowLeft == NULL) ? false : CU::isSameSlice(*cuCurr, *cuBelowLeft);
    isBelowRightAvail = (cuBelowRight == NULL) ? false : CU::isSameSlice(*cuCurr, *cuBelowRight);
  } else {
    isLeftAvail = (cuLeft != NULL);
    isAboveAvail = (cuAbove != NULL);
    isRightAvail = (cuRight != NULL);
    isBelowAvail = (cuBelow != NULL);
    isAboveLeftAvail = (cuAboveLeft != NULL);
    isAboveRightAvail = (cuAboveRight != NULL);
    isBelowLeftAvail = (cuBelowLeft != NULL);
    isBelowRightAvail = (cuBelowRight != NULL);
  }

  // check cross tile flags
  const bool isLoopFilterAcrossTilePPS = cs.pps->getLoopFilterAcrossTilesEnabledFlag();
  if (!isLoopFilterAcrossTilePPS) {
    isLeftAvail = (!isLeftAvail) ? false : CU::isSameTile(*cuCurr, *cuLeft);
    isAboveAvail = (!isAboveAvail) ? false : CU::isSameTile(*cuCurr, *cuAbove);
    isRightAvail = (!isRightAvail) ? false : CU::isSameTile(*cuCurr, *cuRight);
    isBelowAvail = (!isBelowAvail) ? false : CU::isSameTile(*cuCurr, *cuBelow);
    isAboveLeftAvail = (!isAboveLeftAvail) ? false : CU::isSameTile(*cuCurr, *cuAboveLeft);
    isAboveRightAvail = (!isAboveRightAvail) ? false : CU::isSameTile(*cuCurr, *cuAboveRight);
    isBelowLeftAvail = (!isBelowLeftAvail) ? false : CU::isSameTile(*cuCurr, *cuBelowLeft);
    isBelowRightAvail = (!isBelowRightAvail) ? false : CU::isSameTile(*cuCurr, *cuBelowRight);
  }
  // check cross subpic flags
  const SPS& sps = *cs.sps;

  if (sps.getSubPicInfoPresentFlag()) {
    const SubPic& curSubPic = cs.pps->getSubPicFromCU(*cuCurr);
    if (!curSubPic.getloopFilterAcrossSubPicEnabledFlag()) {
      isLeftAvail = (!isLeftAvail) ? false : CU::isSameSubPic(*cuCurr, *cuLeft);
      isAboveAvail = (!isAboveAvail) ? false : CU::isSameSubPic(*cuCurr, *cuAbove);
      isRightAvail = (!isRightAvail) ? false : CU::isSameSubPic(*cuCurr, *cuRight);
      isBelowAvail = (!isBelowAvail) ? false : CU::isSameSubPic(*cuCurr, *cuBelow);
      isAboveLeftAvail = (!isAboveLeftAvail) ? false : CU::isSameSubPic(*cuCurr, *cuAboveLeft);
      isAboveRightAvail = (!isAboveRightAvail) ? false : CU::isSameSubPic(*cuCurr, *cuAboveRight);
      isBelowLeftAvail = (!isBelowLeftAvail) ? false : CU::isSameSubPic(*cuCurr, *cuBelowLeft);
      isBelowRightAvail = (!isBelowRightAvail) ? false : CU::isSameSubPic(*cuCurr, *cuBelowRight);
    }
  }
}

bool SampleAdaptiveOffset::isCrossedByVirtualBoundaries(const PicHeader* picHeader, const Area& area,
                                                        int& numHorVirBndry, int& numVerVirBndry, int horVirBndryPos[],
                                                        int verVirBndryPos[]) {
  numHorVirBndry = 0;
  numVerVirBndry = 0;
  if (!picHeader->getVirtualBoundariesPresentFlag()) {
    return false;
  }

  for (int i = 0; i < picHeader->getNumHorVirtualBoundaries(); i++) {
    if (area.y <= picHeader->getVirtualBoundariesPosY(i) &&
        picHeader->getVirtualBoundariesPosY(i) <= area.y + area.height) {
      horVirBndryPos[numHorVirBndry++] = picHeader->getVirtualBoundariesPosY(i);
    }
  }
  for (int i = 0; i < picHeader->getNumVerVirtualBoundaries(); i++) {
    if (area.x <= picHeader->getVirtualBoundariesPosX(i) &&
        picHeader->getVirtualBoundariesPosX(i) <= area.x + area.width) {
      verVirBndryPos[numVerVirBndry++] = picHeader->getVirtualBoundariesPosX(i);
    }
  }

  return numHorVirBndry > 0 || numVerVirBndry > 0;
}

bool SampleAdaptiveOffset::isProcessDisabled(int xPos, int yPos, int numVerVirBndry, int numHorVirBndry,
                                             int verVirBndryPos[], int horVirBndryPos[]) {
  for (int i = 0; i < numVerVirBndry; i++) {
    if ((xPos == verVirBndryPos[i]) || (xPos == verVirBndryPos[i] - 1)) {
      return true;
    }
  }
  for (int i = 0; i < numHorVirBndry; i++) {
    if ((yPos == horVirBndryPos[i]) || (yPos == horVirBndryPos[i] - 1)) {
      return true;
    }
  }

  return false;
}
//! \}
