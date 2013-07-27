// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_IMAGE_HOLEFILLING_H
#define NV_IMAGE_HOLEFILLING_H

#include <nvimage/nvimage.h>

namespace nv 
{
    class FloatImage;
    class BitMap;

    // The coverage index indicates the channel that should not be extrapolated/filled.

    // These functions leave the bitmap intact:
    NVIMAGE_API void fillVoronoi(FloatImage * img, const BitMap * bmap, int coverageIndex = -1);
    NVIMAGE_API void fillPullPush(FloatImage * img, const BitMap * bmap, int coverageIndex = -1);

    // These other functions update the bitmap in each extrapolation step:
    NVIMAGE_API void fillExtrapolate(int passCount, FloatImage * img, BitMap * bmap, int coverageIndex = -1);
    NVIMAGE_API void fillQuadraticExtrapolate(int passCount, FloatImage * img, BitMap * bmap, int coverageIndex = -1);

} // nv namespace

#endif // NV_IMAGE_HOLEFILLING_H
