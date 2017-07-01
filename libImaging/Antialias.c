/*
 * The Python Imaging Library
 * $Id$
 *
 * pilopen antialiasing support
 *
 * history:
 * 2002-03-09 fl  Created (for PIL 1.1.3)
 * 2002-03-10 fl  Added support for mode "F"
 *
 * Copyright (c) 1997-2002 by Secret Labs AB
 *
 * See the README file for information on usage and redistribution.
 */

#include "Imaging.h"

#include <math.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <smmintrin.h>

#ifdef __AVX2__
    #include "immintrin.h"
#endif


/* resampling filters (from antialias.py) */

struct filter {
    float (*filter)(float x);
    float support;
};

static inline float sinc_filter(float x)
{
    if (x == 0.0)
        return 1.0;
    x = x * M_PI;
    return sin(x) / x;
}

static inline float antialias_filter(float x)
{
    /* lanczos (truncated sinc) */
    if (-3.0 <= x && x < 3.0)
        return sinc_filter(x) * sinc_filter(x/3);
    return 0.0;
}

static struct filter ANTIALIAS = { antialias_filter, 3.0 };

static inline float nearest_filter(float x)
{
    if (-0.5 <= x && x < 0.5)
        return 1.0;
    return 0.0;
}

static struct filter NEAREST = { nearest_filter, 0.5 };

static inline float bilinear_filter(float x)
{
    if (x < 0.0)
        x = -x;
    if (x < 1.0)
        return 1.0-x;
    return 0.0;
}

static struct filter BILINEAR = { bilinear_filter, 1.0 };

static inline float bicubic_filter(float x)
{
    /* FIXME: double-check this algorithm */
    /* FIXME: for best results, "a" should be -0.5 to -1.0, but we'll
       set it to zero for now, to match the 1.1 magnifying filter */
#define a 0.0
    if (x < 0.0)
        x = -x;
    if (x < 1.0)
        return (((a + 2.0) * x) - (a + 3.0)) * x*x + 1;
    if (x < 2.0)
        return (((a * x) - 5*a) * x + 8) * x - 4*a;
    return 0.0;
#undef a
}

static struct filter BICUBIC = { bicubic_filter, 2.0 };


/* 8 bits for result. Filter can have negative areas.
   In one cases the sum of the coefficients will be negative,
   in the other it will be more than 1.0. That is why we need
   two extra bits for overflow and int type. */
#define PRECISION_BITS (32 - 8 - 2)


void
ImagingResampleHorizontalConvolution8u(UINT32 *lineOut, UINT32 *lineIn,
    int xsize, int *xbounds, float *kk, int kmax)
{
    int xmin, xmax, xx, x;
    float *k;

    for (xx = 0; xx < xsize; xx++) {
        xmin = xbounds[xx * 2 + 0];
        xmax = xbounds[xx * 2 + 1];
        k = &kk[xx * kmax];
        x = xmin;
#ifdef __AVX2__
        __m256 sss256 = _mm256_set1_ps(0.25);
        for (; x < xmax - 1; x += 2) {
            __m256 mmk = _mm256_set1_ps(k[x - xmin]);
            mmk = _mm256_insertf128_ps(mmk, _mm_set1_ps(k[x - xmin + 1]), 1);
            __m256i pix = _mm256_cvtepu8_epi32(*(__m128i *) &lineIn[x]);
            __m256 mul = _mm256_mul_ps(_mm256_cvtepi32_ps(pix), mmk);
            sss256 = _mm256_add_ps(sss256, mul);
        }
        __m128 sss = _mm_add_ps(
            _mm256_castps256_ps128(sss256),
            _mm256_extractf128_ps(sss256, 1));
#else
        __m128 sss = _mm_set1_ps(0.5);
#endif
        for (; x < xmax; x++) {
            __m128i pix = _mm_cvtepu8_epi32(*(__m128i *) &lineIn[x]);
            __m128 mmk = _mm_set1_ps(k[x - xmin]);
            __m128 mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);
        }

        __m128i ssi = _mm_cvtps_epi32(sss);
        ssi = _mm_packs_epi32(ssi, ssi);
        lineOut[xx] = _mm_cvtsi128_si32(_mm_packus_epi16(ssi, ssi));
    }
}


void
ImagingResampleVerticalConvolution8u(UINT32 *lineOut, Imaging imIn,
    int ymin, int ymax, float *k)
{
    int y, xx = 0;

#ifdef __AVX2__
    for (; xx < imIn->xsize - 1; xx += 2) {
        __m256 sss = _mm256_set1_ps(0.5);
        for (y = ymin; y < ymax; y++) {
            __m256i pix = _mm256_cvtepu8_epi32(*(__m128i *) &imIn->image32[y][xx]);
            __m256 mmk = _mm256_set1_ps(k[y - ymin]);
            __m256 mul = _mm256_mul_ps(_mm256_cvtepi32_ps(pix), mmk);
            sss = _mm256_add_ps(sss, mul);
        }
        __m256i ssi = _mm256_cvtps_epi32(sss);
        ssi = _mm256_packs_epi32(ssi, ssi);
        ssi = _mm256_packus_epi16(ssi, ssi);
        _mm_storel_epi64((__m128i *) &lineOut[xx], _mm256_castsi256_si128(ssi));
    }
#endif

    for (; xx < imIn->xsize; xx++) {
        __m128 sss = _mm_set1_ps(0.5);
        for (y = ymin; y < ymax; y++) {
            __m128i pix = _mm_cvtepu8_epi32(*(__m128i *) &imIn->image32[y][xx]);
            __m128 mmk = _mm_set1_ps(k[y - ymin]);
            __m128 mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);
        }
        __m128i ssi = _mm_cvtps_epi32(sss);
        ssi = _mm_packs_epi32(ssi, ssi);
        lineOut[xx] = _mm_cvtsi128_si32(_mm_packus_epi16(ssi, ssi));
    }
}



Imaging
ImagingStretch(Imaging imOut, Imaging imIn, int filter)
{
    /* FIXME: this is a quick and straightforward translation from a
       python prototype.  might need some further C-ification... */

    ImagingSectionCookie cookie;
    struct filter *filterp;
    float support, scale, filterscale;
    float center, ww, ss;
    int ymin, ymax, xmin, xmax;
    int kmax, xx, yy, x, y;
    int *intk, *intkk;
    float *k, *kk;
    int *xbounds;

    /* check modes */
    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();

    /* check filter */
    switch (filter) {
    case IMAGING_TRANSFORM_NEAREST:
        filterp = &NEAREST;
        break;
    case IMAGING_TRANSFORM_ANTIALIAS:
        filterp = &ANTIALIAS;
        break;
    case IMAGING_TRANSFORM_BILINEAR:
        filterp = &BILINEAR;
        break;
    case IMAGING_TRANSFORM_BICUBIC:
        filterp = &BICUBIC;
        break;
    default:
        return (Imaging) ImagingError_ValueError(
            "unsupported resampling filter"
            );
    }

    if (imIn->ysize == imOut->ysize) {
        /* prepare for horizontal stretch */
        filterscale = scale = (float) imIn->xsize / imOut->xsize;
    } else if (imIn->xsize == imOut->xsize) {
        /* prepare for vertical stretch */
        filterscale = scale = (float) imIn->ysize / imOut->ysize;
    } else
        return (Imaging) ImagingError_Mismatch();

    /* determine support size (length of resampling filter) */
    support = filterp->support;

    if (filterscale < 1.0) {
        filterscale = 1.0;
    }

    support = support * filterscale;
    kmax = (int) ceil(support) * 2 + 1;

    ImagingSectionEnter(&cookie);
    if (imIn->xsize == imOut->xsize) {
        /* vertical stretch */

        /* coefficient buffer (with rounding safety margin) */
        k = malloc(kmax * sizeof(float));
        if (!k)
            return (Imaging) ImagingError_MemoryError();
        intk = malloc(kmax * sizeof(int));
        if (!intk) {
            free(k);
            return (Imaging) ImagingError_MemoryError();
        }

        for (yy = 0; yy < imOut->ysize; yy++) {
            center = (yy + 0.5) * scale;
            ww = 0.0;
            ss = 1.0 / filterscale;
            /* calculate filter weights */
            ymin = (int) floor(center - support);
            if (ymin < 0)
                ymin = 0;
            ymax = (int) ceil(center + support);
            if (ymax > imIn->ysize)
                ymax = imIn->ysize;
            for (y = ymin; y < ymax; y++) {
                float w = filterp->filter((y - center + 0.5) * ss) * ss;
                k[y - ymin] = w;
                ww = ww + w;
            }
            for (y = 0; y < ymax - ymin; y++) {
                k[y] /= ww;
                intk[y] = (int) (k[y] * (1 << PRECISION_BITS));
            }

            switch(imIn->type) {
            case IMAGING_TYPE_UINT8:
                ImagingResampleVerticalConvolution8u(
                    (UINT32 *) imOut->image32[yy], imIn,
                    ymin, ymax, k
                );
                break;
            }
        }
        free(k);
        free(intk);
    } else {
        kk = malloc(imOut->xsize * kmax * sizeof(float));
        if (!kk)
            return (Imaging) ImagingError_MemoryError();
        intkk = malloc(imOut->xsize * kmax * sizeof(int));
        if (!intkk) {
            free(kk);
            return (Imaging) ImagingError_MemoryError();
        }
        xbounds = malloc(imOut->xsize * 2 * sizeof(int));
        if ( ! xbounds) {
            free(kk);
            free(intkk);
            return 0;
        }

        /* horizontal stretch */
        for (xx = 0; xx < imOut->xsize; xx++) {
            center = (xx + 0.5) * scale;
            ww = 0.0;
            ss = 1.0 / filterscale;
            xmin = (int) floor(center - support);
            if (xmin < 0)
                xmin = 0;
            xmax = (int) ceil(center + support);
            if (xmax > imIn->xsize)
                xmax = imIn->xsize;
            k = &kk[xx * kmax];
            intk = &intkk[xx * kmax];
            for (x = xmin; x < xmax; x++) {
                float w = filterp->filter((x - center + 0.5) * ss) * ss;
                k[x - xmin] = w;
                ww = ww + w;
            }
            for (x = 0; x < xmax - xmin; x++) {
                k[x] /= ww;
                intk[x] = (int) (k[x] * (1 << PRECISION_BITS));
            }
            xbounds[xx * 2 + 0] = xmin;
            xbounds[xx * 2 + 1] = xmax;
        }
        for (yy = 0; yy < imOut->ysize; yy++) {
                switch(imIn->type) {
                case IMAGING_TYPE_UINT8:
                    ImagingResampleHorizontalConvolution8u(
                        (UINT32 *) imOut->image32[yy],
                        (UINT32 *) imIn->image32[yy],
                        imOut->xsize, xbounds, kk, kmax
                    );
                    break;
                }
        }
        free(kk);
        free(intkk);
        free(xbounds);
    }
    ImagingSectionLeave(&cookie);

    return imOut;
}
