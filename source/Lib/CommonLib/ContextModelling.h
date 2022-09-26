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

/** \file     ContextModelling.h
 *  \brief    Classes providing probability descriptions and contexts (header)
 */

#ifndef __CONTEXTMODELLING__
#define __CONTEXTMODELLING__

#include "CommonDef.h"
#include "Contexts.h"
#include "Slice.h"
#include "Unit.h"
#include "UnitPartitioner.h"

#include <bitset>

#include "CommonLib/dtrace_next.h"
class CUCtx {
 public:
  CUCtx() : isDQPCoded(false), isChromaQpAdjCoded(false), qgStart(false), lfnstLastScanPos(false) {
    violatesLfnstConstrained[CHANNEL_TYPE_LUMA] = false;
    violatesLfnstConstrained[CHANNEL_TYPE_CHROMA] = false;
    violatesMtsCoeffConstraint = false;
    mtsLastScanPos = false;
  }
  CUCtx(int _qp) : isDQPCoded(false), isChromaQpAdjCoded(false), qgStart(false), lfnstLastScanPos(false), qp(_qp) {
    violatesLfnstConstrained[CHANNEL_TYPE_LUMA] = false;
    violatesLfnstConstrained[CHANNEL_TYPE_CHROMA] = false;
    violatesMtsCoeffConstraint = false;
    mtsLastScanPos = false;
  }
  ~CUCtx() {}

 public:
  bool isDQPCoded;
  bool isChromaQpAdjCoded;
  bool qgStart;
  bool lfnstLastScanPos;
  int8_t qp;  // used as a previous(last) QP and for QP prediction
  bool violatesLfnstConstrained[MAX_NUM_CHANNEL_TYPE];
  bool violatesMtsCoeffConstraint;
  bool mtsLastScanPos;
};

class MergeCtx {
 public:
  MergeCtx() : numValidMergeCand(0) { memset(mrgTypeNeighbours, 0, sizeof(mrgTypeNeighbours)); }
  ~MergeCtx() {}

 public:
  MvField mvFieldNeighbours[MRG_MAX_NUM_CANDS << 1];  // double length for mv of both lists
  uint8_t BcwIdx[MRG_MAX_NUM_CANDS];
  unsigned char interDirNeighbours[MRG_MAX_NUM_CANDS];
  MergeType mrgTypeNeighbours[MRG_MAX_NUM_CANDS];
  int numValidMergeCand;

  MotionBuf subPuMvpMiBuf;
  MvField mmvdBaseMv[MMVD_BASE_MV_NUM][2];

  void setMmvdMergeCandiInfo(PredictionUnit& pu, int candIdx);
  bool mmvdUseAltHpelIf[MMVD_BASE_MV_NUM];
  bool useAltHpelIf[MRG_MAX_NUM_CANDS];
  void setMergeInfo(PredictionUnit& pu, int candIdx);
  void init() {
    numValidMergeCand = 0;
    memset(mrgTypeNeighbours, 0, sizeof(mrgTypeNeighbours));
  }
};

class AffineMergeCtx {
 public:
  AffineMergeCtx() : numValidMergeCand(0) {
    for (unsigned i = 0; i < AFFINE_MRG_MAX_NUM_CANDS; i++) affineType[i] = AFFINEMODEL_4PARAM;
  }
  ~AffineMergeCtx() {}

 public:
  MvField mvFieldNeighbours[AFFINE_MRG_MAX_NUM_CANDS << 1][3];  // double length for mv of both lists
  unsigned char interDirNeighbours[AFFINE_MRG_MAX_NUM_CANDS];
  AffineModel affineType[AFFINE_MRG_MAX_NUM_CANDS];
  uint8_t BcwIdx[AFFINE_MRG_MAX_NUM_CANDS];
  int numValidMergeCand;
  int maxNumMergeCand;

  MergeCtx* mrgCtx;
  MergeType mergeType[AFFINE_MRG_MAX_NUM_CANDS];
};

namespace DeriveCtx {
void CtxSplit(const CodingStructure& cs, Partitioner& partitioner, unsigned& ctxSpl, unsigned& ctxQt, unsigned& ctxHv,
              unsigned& ctxHorBt, unsigned& ctxVerBt, bool* canSplit = nullptr);
unsigned CtxModeConsFlag(const CodingStructure& cs, Partitioner& partitioner);
unsigned CtxQtCbf(const ComponentID compID, const bool prevCbCbf = false, const int ispIdx = 0);
unsigned CtxInterDir(const PredictionUnit& pu);
unsigned CtxSkipFlag(const CodingUnit& cu);
unsigned CtxAffineFlag(const CodingUnit& cu);
unsigned CtxPredModeFlag(const CodingUnit& cu);
unsigned CtxIBCFlag(const CodingUnit& cu);
unsigned CtxMipFlag(const CodingUnit& cu);
}  // namespace DeriveCtx

#endif  // __CONTEXTMODELLING__
