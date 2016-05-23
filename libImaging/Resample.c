/*
 * The Python Imaging Library
 * $Id$
 *
 * Pillow image resampling support
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

#if defined(__AVX2__)
    #include <immintrin.h>
#endif


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

static inline float lanczos_filter(float x)
{
    /* truncated sinc */
    if (-3.0 <= x && x < 3.0)
        return sinc_filter(x) * sinc_filter(x/3);
    return 0.0;
}

static struct filter LANCZOS = { lanczos_filter, 3.0 };

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
    /* https://en.wikipedia.org/wiki/Bicubic_interpolation#Bicubic_convolution_algorithm */
#define a -0.5
    if (x < 0.0)
        x = -x;
    if (x < 1.0)
        return ((a + 2.0) * x - (a + 3.0)) * x*x + 1;
    if (x < 2.0)
        return (((x - 5) * x + 8) * x - 4) * a;
    return 0.0;
#undef a
}

static struct filter BICUBIC = { bicubic_filter, 2.0 };


static inline UINT8 clip8(float in)
{
    int out = (int) in;
    if (out >= 255)
       return 255;
    if (out <= 0)
        return 0;
    return (UINT8) out;
}


void
ImagingResampleHorizontalConvolution8u(UINT32 *lineOut, UINT32 *lineIn,
    int xsize, int *xbounds, float *kk, int kmax)
{
    int xmin, xmax, xx, x;
    float *k;

    __m128i mmmax = _mm_set1_epi32(255);
    __m128i mmmin = _mm_set1_epi32(0);
    __m128i shiftmask = _mm_set_epi8(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,12,8,4,0);

    for (xx = 0; xx < xsize; xx++) {
        __m128i pix;
        __m128 mmk, mul, sss;
        xmin = xbounds[xx * 2 + 0];
        xmax = xbounds[xx * 2 + 1];
        k = &kk[xx * kmax];
        x = 0;

#if defined(__AVX2__)

        __m256i pix256;
        __m256 mmk256, mul256;
        __m256 sss256 = _mm256_setzero_ps();

        for (; x < xmax - 7; x += 8) {
            __m256i source = _mm256_loadu_si256((__m256i *) &lineIn[x + xmin]);
            __m256 sourcek = _mm256_loadu_ps(&k[x]);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0,
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            mmk256 = _mm256_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(0, 0, 0, 0));
            mul256 = _mm256_mul_ps(_mm256_cvtepi32_ps(pix256), mmk256);
            sss256 = _mm256_add_ps(sss256, mul256);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4,
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            mmk256 = _mm256_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(1, 1, 1, 1));
            mul256 = _mm256_mul_ps(_mm256_cvtepi32_ps(pix256), mmk256);
            sss256 = _mm256_add_ps(sss256, mul256);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8,
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            mmk256 = _mm256_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(2, 2, 2, 2));
            mul256 = _mm256_mul_ps(_mm256_cvtepi32_ps(pix256), mmk256);
            sss256 = _mm256_add_ps(sss256, mul256);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12,
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            mmk256 = _mm256_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(3, 3, 3, 3));
            mul256 = _mm256_mul_ps(_mm256_cvtepi32_ps(pix256), mmk256);
            sss256 = _mm256_add_ps(sss256, mul256);
        }

        for (; x < xmax - 1; x += 2) {
            pix256 = _mm256_cvtepu8_epi32(_mm_set_epi64x(0,  *((uint64_t *)&lineIn[x + xmin]) ));
            mmk256 = _mm256_insertf128_ps(_mm256_set1_ps(k[x]), _mm_set1_ps(k[x + 1]), 1);
            mul256 = _mm256_mul_ps(_mm256_cvtepi32_ps(pix256), mmk256);
            sss256 = _mm256_add_ps(sss256, mul256);
        }

        sss = _mm_add_ps(_mm256_extractf128_ps(sss256, 0),
            _mm256_extractf128_ps(sss256, 1));

#else

        sss = _mm_setzero_ps();
        for (; x < xmax - 3; x += 4) {
            __m128i source = _mm_loadu_si128((__m128i *) &lineIn[x + xmin]);
            __m128 sourcek = _mm_loadu_ps(&k[x]);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            mmk = _mm_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(0, 0, 0, 0));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            mmk = _mm_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(1, 1, 1, 1));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            mmk = _mm_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(2, 2, 2, 2));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            mmk = _mm_shuffle_ps(sourcek, sourcek, _MM_SHUFFLE(3, 3, 3, 3));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);
        }

#endif

        for (; x < xmax; x++) {
            pix = _mm_cvtepu8_epi32(*(__m128i *) &lineIn[x + xmin]);
            mmk = _mm_set1_ps(k[x]);
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);
        }

        __m128i ssi = _mm_cvtps_epi32(sss);
        ssi = _mm_max_epi32(mmmin, _mm_min_epi32(mmmax, ssi));
        lineOut[xx] = _mm_cvtsi128_si32(_mm_shuffle_epi8(ssi, shiftmask));
    }
}

void
ImagingResampleVerticalConvolution8u(UINT32 *lineOut, Imaging imIn,
    int xmin, int xmax, float *k)
{
    int x;
    int xx = 0;
    int xsize = imIn->xsize;

    __m128i mmmax = _mm_set1_epi32(255);
    __m128i mmmin = _mm_set1_epi32(0);
    __m128i shiftmask = _mm_set_epi8(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,12,8,4,0);
    __m128i pix;
    __m128 mmk, mul;

    for (; xx < xsize - 3; xx += 4) {
        __m128i ssi0, ssi1;
        __m128 sss0 = _mm_setzero_ps();
        __m128 sss1 = _mm_setzero_ps();
        __m128 sss2 = _mm_setzero_ps();
        __m128 sss3 = _mm_setzero_ps();

        for (x = 0; x < xmax; x++) {
            __m128i source = _mm_loadu_si128((__m128i *) &imIn->image32[x + xmin][xx]);
            mmk = _mm_set1_ps(k[x]);
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss0 = _mm_add_ps(sss0, mul);
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss1 = _mm_add_ps(sss1, mul);
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss2 = _mm_add_ps(sss2, mul);
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(-1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss3 = _mm_add_ps(sss3, mul);
        }
        ssi0 = _mm_packus_epi32(_mm_cvtps_epi32(sss0), _mm_cvtps_epi32(sss1));
        ssi1 = _mm_packus_epi32(_mm_cvtps_epi32(sss2), _mm_cvtps_epi32(sss3));
        ssi0 = _mm_packus_epi16(ssi0, ssi1);
        _mm_storeu_si128((__m128i *) &lineOut[xx], ssi0);
    }

    for (; xx < xsize; xx++) {
        __m128i ssi;
        __m128 sss = _mm_setzero_ps();
        for (x = 0; x < xmax; x++) {
            pix = _mm_cvtepu8_epi32(*(__m128i *) &imIn->image32[x + xmin][xx]);
            mmk = _mm_set1_ps(k[x]);
            mul = _mm_mul_ps(_mm_cvtepi32_ps(pix), mmk);
            sss = _mm_add_ps(sss, mul);
        }
        ssi = _mm_cvtps_epi32(sss);
        ssi = _mm_max_epi32(mmmin, _mm_min_epi32(mmmax, ssi));
        lineOut[xx] = _mm_cvtsi128_si32(_mm_shuffle_epi8(ssi, shiftmask));
    }
}


/* This is work around bug in GCC prior 4.9 in 64-bit mode.
   GCC generates code with partial dependency which 3 times slower.
   See: http://stackoverflow.com/a/26588074/253146 */
#if defined(__x86_64__) && defined(__SSE__) &&  ! defined(__NO_INLINE__) && \
    ! defined(__clang__) && defined(GCC_VERSION) && (GCC_VERSION < 40900)
static float __attribute__((always_inline)) i2f(int v) {
    float x;
    __asm__("xorps %0, %0; cvtsi2ss %1, %0" : "=X"(x) : "r"(v) );
    return x;
}
#else
static float inline i2f(int v) { return (float) v; }
#endif


Imaging
ImagingResampleHorizontal(Imaging imIn, int xsize, int filter)
{
    ImagingSectionCookie cookie;
    Imaging imOut;
    struct filter *filterp;
    float support, scale, filterscale;
    float center, ww, ss;
    int xx, yy, x, kmax, xmin, xmax;
    int *xbounds;
    float *k, *kk;

    /* check filter */
    switch (filter) {
    case IMAGING_TRANSFORM_LANCZOS:
        filterp = &LANCZOS;
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

    /* prepare for horizontal stretch */
    filterscale = scale = (float) imIn->xsize / xsize;

    /* determine support size (length of resampling filter) */
    support = filterp->support;

    if (filterscale < 1.0) {
        filterscale = 1.0;
    }

    support = support * filterscale;

    /* maximum number of coofs */
    kmax = (int) ceil(support) * 2 + 1;

    // check for overflow
    if (kmax > 0 && xsize > SIZE_MAX / kmax)
        return (Imaging) ImagingError_MemoryError();

    // sizeof(float) should be greater than 0
    if (xsize * kmax > SIZE_MAX / sizeof(float))
        return (Imaging) ImagingError_MemoryError();

    /* coefficient buffer */
    kk = malloc(xsize * kmax * sizeof(float));
    if ( ! kk)
        return (Imaging) ImagingError_MemoryError();

    // sizeof(int) should be greater than 0 as well
    if (xsize > SIZE_MAX / (2 * sizeof(int)))
        return (Imaging) ImagingError_MemoryError();

    xbounds = malloc(xsize * 2 * sizeof(int));
    if ( ! xbounds) {
        free(kk);
        return (Imaging) ImagingError_MemoryError();
    }

    for (xx = 0; xx < xsize; xx++) {
        k = &kk[xx * kmax];
        center = (xx + 0.5) * scale;
        ww = 0.0;
        ss = 1.0 / filterscale;
        xmin = (int) floor(center - support);
        if (xmin < 0)
            xmin = 0;
        xmax = (int) ceil(center + support);
        if (xmax > imIn->xsize)
            xmax = imIn->xsize;
        for (x = xmin; x < xmax; x++) {
            float w = filterp->filter((x - center + 0.5) * ss) * ss;
            k[x - xmin] = w;
            ww += w;
        }
        for (x = 0; x < xmax - xmin; x++) {
            if (ww != 0.0)
                k[x] /= ww;
        }
        xbounds[xx * 2 + 0] = xmin;
        xbounds[xx * 2 + 1] = xmax - xmin;
    }

    imOut = ImagingNew(imIn->mode, xsize, imIn->ysize);
    if ( ! imOut) {
        free(kk);
        free(xbounds);
        return NULL;
    }

    ImagingSectionEnter(&cookie);
    /* horizontal stretch */
    for (yy = 0; yy < imOut->ysize; yy++) {
        if (imIn->image8) {
            /* 8-bit grayscale */
            for (xx = 0; xx < xsize; xx++) {
                xmin = xbounds[xx * 2 + 0];
                xmax = xbounds[xx * 2 + 1];
                k = &kk[xx * kmax];
                ss = 0.5;
                for (x = 0; x < xmax; x++)
                    ss += i2f(imIn->image8[yy][x + xmin]) * k[x];
                imOut->image8[yy][xx] = clip8(ss);
            }
        } else {
            switch(imIn->type) {
            case IMAGING_TYPE_UINT8:
                /* n-bit grayscale */
                ImagingResampleHorizontalConvolution8u(
                    (UINT32 *) imOut->image32[yy],
                    (UINT32 *) imIn->image32[yy],
                    xsize, xbounds, kk, kmax
                );
                break;
            case IMAGING_TYPE_INT32:
                /* 32-bit integer */
                for (xx = 0; xx < xsize; xx++) {
                    xmin = xbounds[xx * 2 + 0];
                    xmax = xbounds[xx * 2 + 1];
                    k = &kk[xx * kmax];
                    ss = 0.0;
                    for (x = 0; x < xmax; x++)
                        ss += i2f(IMAGING_PIXEL_I(imIn, x + xmin, yy)) * k[x];
                    IMAGING_PIXEL_I(imOut, xx, yy) = (int) ss;
                }
                break;
            case IMAGING_TYPE_FLOAT32:
                /* 32-bit float */
                for (xx = 0; xx < xsize; xx++) {
                    xmin = xbounds[xx * 2 + 0];
                    xmax = xbounds[xx * 2 + 1];
                    k = &kk[xx * kmax];
                    ss = 0.0;
                    for (x = 0; x < xmax; x++)
                        ss += IMAGING_PIXEL_F(imIn, x + xmin, yy) * k[x];
                    IMAGING_PIXEL_F(imOut, xx, yy) = ss;
                }
                break;
            }
        }
    }
    ImagingSectionLeave(&cookie);
    free(kk);
    free(xbounds);
    return imOut;
}


Imaging
ImagingResampleVertical(Imaging imIn, int ysize, int filter)
{
    ImagingSectionCookie cookie;
    Imaging imOut;
    struct filter *filterp;
    float support, scale, filterscale;
    float center, ww, ss;
    int xx, yy, x, kmax, xmin, xmax;
    float *k;

    /* check filter */
    switch (filter) {
    case IMAGING_TRANSFORM_LANCZOS:
        filterp = &LANCZOS;
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

    /* prepare for horizontal stretch */
    filterscale = scale = (float) imIn->ysize / ysize;

    /* determine support size (length of resampling filter) */
    support = filterp->support;

    if (filterscale < 1.0) {
        filterscale = 1.0;
    }

    support = support * filterscale;

    /* maximum number of coofs */
    kmax = (int) ceil(support) * 2 + 1;

    // sizeof(float) should be greater than 0
    if (kmax > SIZE_MAX / sizeof(float))
        return (Imaging) ImagingError_MemoryError();

    /* coefficient buffer (with rounding safety margin) */
    k = malloc(kmax * sizeof(float));
    if ( ! k)
        return (Imaging) ImagingError_MemoryError();

    imOut = ImagingNew(imIn->mode, imIn->xsize, ysize);
    if ( ! imOut) {
        free(k);
        return NULL;
    }

    ImagingSectionEnter(&cookie);
    /* horizontal stretch */
    for (yy = 0; yy < ysize; yy++) {
        center = (yy + 0.5) * scale;
        ww = 0.0;
        ss = 1.0 / filterscale;
        xmin = (int) floor(center - support);
        if (xmin < 0)
            xmin = 0;
        xmax = (int) ceil(center + support);
        if (xmax > imIn->ysize)
            xmax = imIn->ysize;
        for (x = xmin; x < xmax; x++) {
            float w = filterp->filter((x - center + 0.5) * ss) * ss;
            k[x - xmin] = w;
            ww += w;
        }
        for (x = 0; x < xmax - xmin; x++) {
            if (ww != 0.0)
                k[x] /= ww;
        }

        if (imIn->image8) {
            /* 8-bit grayscale */
            for (xx = 0; xx < imIn->xsize; xx++) {
                ss = 0.5;
                for (x = xmin; x < xmax; x++)
                    ss += i2f(imIn->image8[x][xx]) * k[x - xmin];
                imOut->image8[yy][xx] = clip8(ss);
            }
        } else {
            switch(imIn->type) {
            case IMAGING_TYPE_UINT8:
                /* n-bit grayscale */
                ImagingResampleVerticalConvolution8u(
                    (UINT32 *) imOut->image32[yy], imIn,
                    xmin, xmax - xmin, k
                );
                break;
            case IMAGING_TYPE_INT32:
                /* 32-bit integer */
                for (xx = 0; xx < imIn->xsize; xx++) {
                    ss = 0.0;
                    for (x = xmin; x < xmax; x++)
                        ss += i2f(IMAGING_PIXEL_I(imIn, xx, x)) * k[x - xmin];
                    IMAGING_PIXEL_I(imOut, xx, yy) = (int) ss;
                }
                break;
            case IMAGING_TYPE_FLOAT32:
                /* 32-bit float */
                for (xx = 0; xx < imIn->xsize; xx++) {
                    ss = 0.0;
                    for (x = xmin; x < xmax; x++)
                        ss += IMAGING_PIXEL_F(imIn, xx, x) * k[x - xmin];
                    IMAGING_PIXEL_F(imOut, xx, yy) = ss;
                }
                break;
            default:
                ImagingSectionLeave(&cookie);
                return (Imaging) ImagingError_ModeError();
            }
        }
    }
    ImagingSectionLeave(&cookie);
    free(k);
    return imOut;
}


Imaging
ImagingResample(Imaging imIn, int xsize, int ysize, int filter)
{
    Imaging imTemp;
    Imaging imOut;

    if (strcmp(imIn->mode, "P") == 0 || strcmp(imIn->mode, "1") == 0)
        return (Imaging) ImagingError_ModeError();

    if (imIn->type == IMAGING_TYPE_SPECIAL)
        return (Imaging) ImagingError_ModeError();

    imTemp = ImagingResampleHorizontal(imIn, xsize, filter);
    if ( ! imTemp)
        return NULL;

    imOut = ImagingResampleVertical(imTemp, ysize, filter);
    ImagingDelete(imTemp);

    return imOut;
}
