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


static inline UINT8 clip8(float in)
{
    int out = (int) in;
    if (out >= 255)
       return 255;
    if (out <= 0)
        return 0;
    return (UINT8) out;
}

#ifdef __GNUC__
    #define GCC_VERSION (__GNUC__ * 10000 \
                       + __GNUC_MINOR__ * 100 \
                       + __GNUC_PATCHLEVEL__)
#endif

/* This is work around bug in GCC prior 4.9 in 64 bit mode.
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
ImagingStretch(Imaging imOut, Imaging imIn, int filter)
{
    /* FIXME: this is a quick and straightforward translation from a
       python prototype.  might need some further C-ification... */

    ImagingSectionCookie cookie;
    struct filter *filterp;
    float support, scale, filterscale;
    float center, ww;
    float ss, ss0, ss1, ss2;
    int ymin, ymax, xmin, xmax;
    int kmax, xx, yy, x, y;
    float *k;
    float *kk;
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
            for (y = ymin; y < ymax; y++) {
                k[y - ymin] /= ww;
            }

            switch(imIn->type) {
            case IMAGING_TYPE_UINT8:
                if (imIn->bands == 4) {
                    // Тело для 4-х каналов
                } else if (imIn->bands == 3) {
                    for (xx = 0; xx < imOut->xsize; xx++) {
                        ss0 = 0.5;
                        ss1 = 0.5;
                        ss2 = 0.5;
                        for (y = ymin; y < ymax; y++) {
                            ss0 = ss0 + i2f((UINT8) imIn->image[y][xx*4+0]) * k[y-ymin];
                            ss1 = ss1 + i2f((UINT8) imIn->image[y][xx*4+1]) * k[y-ymin];
                            ss2 = ss2 + i2f((UINT8) imIn->image[y][xx*4+2]) * k[y-ymin];
                        }
                        imOut->image[yy][xx*4+0] = clip8(ss0);
                        imOut->image[yy][xx*4+1] = clip8(ss1);
                        imOut->image[yy][xx*4+2] = clip8(ss2);
                    }
                }
                break;
            }
        }
        free(k);
    } else {
        kk = malloc(imOut->xsize * kmax * sizeof(float));
        if (!kk)
            return (Imaging) ImagingError_MemoryError();
        xbounds = malloc(imOut->xsize * 2 * sizeof(int));
        if ( ! xbounds) {
            free(kk);
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
            for (x = xmin; x < xmax; x++) {
                float w = filterp->filter((x - center + 0.5) * ss) * ss;
                k[x - xmin] = w;
                ww = ww + w;
            }
            for (x = xmin; x < xmax; x++) {
                k[x - xmin] /= ww;
            }
            xbounds[xx * 2 + 0] = xmin;
            xbounds[xx * 2 + 1] = xmax;
        }
        for (yy = 0; yy < imOut->ysize; yy++) {
                switch(imIn->type) {
                case IMAGING_TYPE_UINT8:
                    if (imIn->bands == 4) {
                        // Тело для 4-х каналов
                    } else if (imIn->bands == 3) {
                        for (xx = 0; xx < imOut->xsize; xx++) {
                            k = &kk[xx * kmax];
                            xmin = xbounds[xx * 2 + 0];
                            xmax = xbounds[xx * 2 + 1];
                            ss0 = 0.5;
                            ss1 = 0.5;
                            ss2 = 0.5;
                            for (x = xmin; x < xmax; x++){
                                ss0 = ss0 + i2f((UINT8) imIn->image[yy][x*4+0]) * k[x - xmin];
                                ss1 = ss1 + i2f((UINT8) imIn->image[yy][x*4+1]) * k[x - xmin];
                                ss2 = ss2 + i2f((UINT8) imIn->image[yy][x*4+2]) * k[x - xmin];
                            }
                            imOut->image[yy][xx*4+0] = clip8(ss0);
                            imOut->image[yy][xx*4+1] = clip8(ss1);
                            imOut->image[yy][xx*4+2] = clip8(ss2);
                        }
                    } else {
                        // Тело для двух и одного канала
                    }
                    break;
                }
        }
        free(kk);
        free(xbounds);
    }
    ImagingSectionLeave(&cookie);

    return imOut;
}
