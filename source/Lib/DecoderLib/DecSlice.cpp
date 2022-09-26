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

/** \file     DecSlice.cpp
    \brief    slice decoder class
*/

#include "DecSlice.h"

#include "CommonLib/SampleAdaptiveOffset.h"
#include "CommonLib/TrQuant.h"
#include "CommonLib/UnitTools.h"
#include "CommonLib/dtrace_next.h"
#include "CommonLib/TimeProfiler.h"
#include "CommonLib/AdaptiveLoopFilter.h"
#include "DecLibRecon.h"
#include "DecLib.h"
#include <vector>

//! \ingroup DecoderLib
//! \{

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void DecSlice::parseSlice(Slice* slice, InputBitstream* bitstream, int threadId, DecLibRecon* recon) {
  PROFILER_SCOPE_AND_STAGE_EXT(1, g_timeProfiler, P_CONTROL_PARSE_DERIVE_LL, *slice->getPic()->cs, CH_L);
  const unsigned numSubstreams = slice->getNumberOfSubstreamSizes() + 1;
  // Table of extracted substreams.
  std::vector<std::shared_ptr<InputBitstream>> ppcSubstreams(numSubstreams);
  for (unsigned idx = 0; idx < numSubstreams; idx++) {
    ppcSubstreams[idx].reset(bitstream->extractSubstream(idx + 1 < numSubstreams ? (slice->getSubstreamSize(idx) << 3)
                                                                                 : bitstream->getNumBitsLeft()));
  }

  const SPS* sps = slice->getSPS();
  Picture* pic = slice->getPic();

  // setup coding structure
  CodingStructure& cs = *pic->cs;
  cs.chromaQpAdj = 0;

  const int startCtuTsAddr = slice->getFirstCtuRsAddrInSlice();
  const unsigned widthInCtus = cs.pcv->widthInCtus;
  const bool wavefrontsEnabled = cs.sps->getEntropyCodingSyncEnabledFlag();
  const bool entryPointPresent = cs.sps->getEntryPointsPresentFlag();

  if (startCtuTsAddr == 0) {
    cs.picture->resizeCtuAlfData(cs.pcv->sizeInCtus);
  }
  AdaptiveLoopFilter::reconstructCoeffAPSs(*slice);

  CABACDecoder cabacDecoder;
  CABACReader& cabacReader = *cabacDecoder.getCABACReader();
  cabacReader.initBitstream(ppcSubstreams[0].get());
  cabacReader.initCtxModels(*slice);

  // if first slice, finish init of the coding structure
  if (startCtuTsAddr == 0) {
    cs.initStructData();
  }

  const uint32_t log2SaoOffsetScaleLuma =
      (uint32_t)std::max(0, cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) - MAX_SAO_TRUNCATED_BITDEPTH);
  const uint32_t log2SaoOffsetScaleChroma =
      (uint32_t)std::max(0, cs.sps->getBitDepth(CHANNEL_TYPE_CHROMA) - MAX_SAO_TRUNCATED_BITDEPTH);
  cabacReader.setUpSAOBitDepth(log2SaoOffsetScaleLuma, log2SaoOffsetScaleChroma);

  // Quantization parameter
  pic->m_prevQP[0] = pic->m_prevQP[1] = slice->getSliceQp();

  CHECK(pic->m_prevQP[0] == std::numeric_limits<int>::max(), "parse slice: Invalid previous QP");

  DTRACE(g_trace_ctx, D_HEADER, "=========== POC: %d ===========\n", slice->getPOC());

  unsigned subStrmId = 0;
  for (unsigned ctuIdx = 0; ctuIdx < slice->getNumCtuInSlice(); ctuIdx++) {
    const unsigned ctuRsAddr = slice->getCtuAddrInSlice(ctuIdx);
    const unsigned ctuXPosInCtus = ctuRsAddr % widthInCtus;
    const unsigned ctuYPosInCtus = ctuRsAddr / widthInCtus;
    const unsigned tileColIdx = slice->getPPS()->ctuToTileCol(ctuXPosInCtus);
    const unsigned tileRowIdx = slice->getPPS()->ctuToTileRow(ctuYPosInCtus);
    const unsigned tileXPosInCtus = slice->getPPS()->getTileColumnBd(tileColIdx);
    const unsigned tileYPosInCtus = slice->getPPS()->getTileRowBd(tileRowIdx);
    const unsigned tileColWidth = slice->getPPS()->getTileColumnWidth(tileColIdx);
    const unsigned tileRowHeight = slice->getPPS()->getTileRowHeight(tileRowIdx);
    const unsigned tileIdx = slice->getPPS()->getTileIdx(ctuXPosInCtus, ctuYPosInCtus);
    const unsigned maxCUSize = sps->getMaxCUWidth();
    Position pos(ctuXPosInCtus * maxCUSize, ctuYPosInCtus * maxCUSize);
    UnitArea ctuArea(cs.area.chromaFormat, Area(pos.x, pos.y, maxCUSize, maxCUSize));
    const SubPic& curSubPic = slice->getPPS()->getSubPicFromPos(pos);
    // padding/restore at slice level
    if (slice->getPPS()->getNumSubPics() >= 2 && curSubPic.getTreatedAsPicFlag() && ctuIdx == 0) {
      int subPicX = (int)curSubPic.getSubPicLeft();
      int subPicY = (int)curSubPic.getSubPicTop();
      int subPicWidth = (int)curSubPic.getSubPicWidthInLumaSample();
      int subPicHeight = (int)curSubPic.getSubPicHeightInLumaSample();
      for (int rlist = REF_PIC_LIST_0; rlist < NUM_REF_PIC_LIST_01; rlist++) {
        int n = slice->getNumRefIdx((RefPicList)rlist);
        for (int idx = 0; idx < n; idx++) {
          Picture* refPic = slice->getRefPic((RefPicList)rlist, idx);
          if (!refPic->getSubPicSaved() && refPic->subPictures.size() > 1) {
            refPic->saveSubPicBorder(refPic->getPOC(), subPicX, subPicY, subPicWidth, subPicHeight);
            refPic->extendSubPicBorder(refPic->getPOC(), subPicX, subPicY, subPicWidth, subPicHeight);
            refPic->setSubPicSaved(true);
          }
        }
      }
    }

    DTRACE_UPDATE(g_trace_ctx, std::make_pair("ctu", ctuRsAddr));

    cabacReader.initBitstream(ppcSubstreams[subStrmId].get());

    // set up CABAC contexts' state for this CTU
    if (ctuXPosInCtus == tileXPosInCtus && ctuYPosInCtus == tileYPosInCtus) {
      if (ctuIdx != 0)  // if it is the first CTU, then the entropy coder has already been reset
      {
        cabacReader.initCtxModels(*slice);
      }
      pic->m_prevQP[0] = pic->m_prevQP[1] = slice->getSliceQp();
    } else if (ctuXPosInCtus == tileXPosInCtus && wavefrontsEnabled) {
      // Synchronize cabac probabilities with upper-right CTU if it's available and at the start of a line.
      if (ctuIdx != 0)  // if it is the first CTU, then the entropy coder has already been reset
      {
        cabacReader.initCtxModels(*slice);
      }
      if (cs.getCURestricted(pos.offset(0, -1), pos, slice->getIndependentSliceIdx(), tileIdx, CH_L)) {
        cabacReader.getCtx() = m_entropyCodingSyncContextState;  // TODO: sync this
      }
      pic->m_prevQP[0] = pic->m_prevQP[1] = slice->getSliceQp();
    }

    // memset( cs.getCtuData( ctuRsAddr ).cuPtr, 0, sizeof( CtuData::cuPtr ) );

    cabacReader.coding_tree_unit(cs, slice, ctuArea, pic->m_prevQP, ctuRsAddr);

    if (recon != nullptr) {
      recon->deriveMV(ctuXPosInCtus, ctuYPosInCtus, threadId);
    }
    if (ctuXPosInCtus == tileXPosInCtus && wavefrontsEnabled) {
      m_entropyCodingSyncContextState = cabacReader.getCtx();  // TODO: sync this
    }

    if (ctuIdx == slice->getNumCtuInSlice() - 1) {
      unsigned binVal = cabacReader.terminating_bit();
      CHECK(!binVal, "parse slice: Expecting a terminating bit");
#if DECODER_CHECK_SUBSTREAM_AND_SLICE_TRAILING_BYTES
      cabacReader.remaining_bytes(false);
#endif
    } else if ((ctuXPosInCtus + 1 == tileXPosInCtus + tileColWidth) &&
               (ctuYPosInCtus + 1 == tileYPosInCtus + tileRowHeight || wavefrontsEnabled)) {
      // The sub-stream/stream should be terminated after this CTU.
      // (end of slice-segment, end of tile, end of wavefront-CTU-row)
      unsigned binVal = cabacReader.terminating_bit();
      CHECK(!binVal, "parse slice: Expecting a terminating bit 2");

      if (entryPointPresent) {
#if DECODER_CHECK_SUBSTREAM_AND_SLICE_TRAILING_BYTES
        cabacReader.remaining_bytes(true);
#endif
        subStrmId++;
      }
    }

    if (slice->getPPS()->getNumSubPics() >= 2 && curSubPic.getTreatedAsPicFlag() &&
        ctuIdx == (slice->getNumCtuInSlice() - 1))
    // for last Ctu in the slice
    {
      int subPicX = (int)curSubPic.getSubPicLeft();
      int subPicY = (int)curSubPic.getSubPicTop();
      int subPicWidth = (int)curSubPic.getSubPicWidthInLumaSample();
      int subPicHeight = (int)curSubPic.getSubPicHeightInLumaSample();
      for (int rlist = REF_PIC_LIST_0; rlist < NUM_REF_PIC_LIST_01; rlist++) {
        int n = slice->getNumRefIdx((RefPicList)rlist);
        for (int idx = 0; idx < n; idx++) {
          Picture* refPic = slice->getRefPic((RefPicList)rlist, idx);
          if (refPic->getSubPicSaved()) {
            refPic->restoreSubPicBorder(refPic->getPOC(), subPicX, subPicY, subPicWidth, subPicHeight);
            refPic->setSubPicSaved(false);
          }
        }
      }
    }
  }
}

//! \}
