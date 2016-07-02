#include "Imaging.h"

#include <math.h>
#include <emmintrin.h>
#include <mmintrin.h>
#include <smmintrin.h>

#if defined(__AVX2__)
    #include <immintrin.h>
#endif


#define ROUND_UP(f) ((int) ((f) >= 0.0 ? (f) + 0.5F : (f) - 0.5F))


struct filter {
    double (*filter)(double x);
    double support;
};

static inline double sinc_filter(double x)
{
    if (x == 0.0)
        return 1.0;
    x = x * M_PI;
    return sin(x) / x;
}

static inline double lanczos_filter(double x)
{
    /* truncated sinc */
    if (-3.0 <= x && x < 3.0)
        return sinc_filter(x) * sinc_filter(x/3);
    return 0.0;
}

static inline double bilinear_filter(double x)
{
    if (x < 0.0)
        x = -x;
    if (x < 1.0)
        return 1.0-x;
    return 0.0;
}

static inline double bicubic_filter(double x)
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

static struct filter LANCZOS = { lanczos_filter, 3.0 };
static struct filter BILINEAR = { bilinear_filter, 1.0 };
static struct filter BICUBIC = { bicubic_filter, 2.0 };



/* 8 bits for result. Filter can have negative areas.
   In one cases the sum of the coefficients will be negative,
   in the other it will be more than 1.0. That is why we need
   two extra bits for overflow and int type. */
#define PRECISION_BITS (32 - 8 - 2)


static inline UINT8 clip8(int in)
{
    if (in >= (1 << PRECISION_BITS << 8))
       return 255;
    if (in <= 0)
        return 0;
    return (UINT8) (in >> PRECISION_BITS);
}


int
ImagingPrecompute(int inSize, int outSize, struct filter *filterp,
                  int **xboundsp, double **kkp) {
    double support, scale, filterscale;
    double center, ww, ss;
    int xx, x, kmax, xmin, xmax;
    int *xbounds;
    double *kk, *k;

    /* prepare for horizontal stretch */
    filterscale = scale = (double) inSize / outSize;
    if (filterscale < 1.0) {
        filterscale = 1.0;
    }

    /* determine support size (length of resampling filter) */
    support = filterp->support * filterscale;

    /* maximum number of coeffs */
    kmax = (int) ceil(support) * 2 + 1;

    // check for overflow
    if (outSize > INT_MAX / (kmax * sizeof(double)))
        return 0;

    /* coefficient buffer */
    /* malloc check ok, overflow checked above */
    kk = malloc(outSize * kmax * sizeof(double));
    if ( ! kk)
        return 0;

    /* malloc check ok, kmax*sizeof(double) > 2*sizeof(int) */
    xbounds = malloc(outSize * 2 * sizeof(int));
    if ( ! xbounds) {
        free(kk);
        return 0;
    }

    for (xx = 0; xx < outSize; xx++) {
        center = (xx + 0.5) * scale;
        ww = 0.0;
        ss = 1.0 / filterscale;
        xmin = (int) floor(center - support);
        if (xmin < 0)
            xmin = 0;
        xmax = (int) ceil(center + support);
        if (xmax > inSize)
            xmax = inSize;
        xmax -= xmin;
        k = &kk[xx * kmax];
        for (x = 0; x < xmax; x++) {
            double w = filterp->filter((x + xmin - center + 0.5) * ss);
            k[x] = w;
            ww += w;
        }
        for (x = 0; x < xmax; x++) {
            if (ww != 0.0)
                k[x] /= ww;
        }
        // Remaining values should stay empty if they are used despite of xmax.
        for (; x < kmax; x++) {
            k[x] = 0;
        }
        xbounds[xx * 2 + 0] = xmin;
        xbounds[xx * 2 + 1] = xmax;
    }
    *xboundsp = xbounds;
    *kkp = kk;
    return kmax;
}


void
ImagingResampleHorizontalConvolution8u(UINT32 *lineOut, UINT32 *lineIn,
    int xsize, int *xbounds, int *kk, int kmax)
{
    int xmin, xmax, xx, x;
    int *k;

    for (xx = 0; xx < xsize; xx++) {
        __m128i sss;
        xmin = xbounds[xx * 2 + 0];
        xmax = xbounds[xx * 2 + 1];
        k = &kk[xx * kmax];
        x = 0;

#if defined(__AVX2__)

        __m256i sss256 = _mm256_set1_epi32(1 << (PRECISION_BITS -2));
        for (; x < xmax - 7; x += 8) {
            __m256i pix256, mmk, mul;
            __m256i source = _mm256_loadu_si256((__m256i *) &lineIn[x + xmin]);
            __m256i sourcek = _mm256_loadu_si256((__m256i *) &k[x]);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0,
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            mmk = _mm256_shuffle_epi32(sourcek, _MM_SHUFFLE(0, 0, 0, 0));
            mul = _mm256_mullo_epi32(pix256, mmk);
            sss256 = _mm256_add_epi32(sss256, mul);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4,
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            mmk = _mm256_shuffle_epi32(sourcek, _MM_SHUFFLE(1, 1, 1, 1));
            mul = _mm256_mullo_epi32(pix256, mmk);
            sss256 = _mm256_add_epi32(sss256, mul);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8,
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            mmk = _mm256_shuffle_epi32(sourcek, _MM_SHUFFLE(2, 2, 2, 2));
            mul = _mm256_mullo_epi32(pix256, mmk);
            sss256 = _mm256_add_epi32(sss256, mul);

            pix256 = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12,
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            mmk = _mm256_shuffle_epi32(sourcek, _MM_SHUFFLE(3, 3, 3, 3));
            mul = _mm256_mullo_epi32(pix256, mmk);
            sss256 = _mm256_add_epi32(sss256, mul);
        }
        for (; x < xmax - 1; x += 2) {
            __m256i mmk = _mm256_set1_epi32(k[x]);
            mmk = _mm256_insertf128_si256(mmk, _mm_set1_epi32(k[x + 1]), 1);
            __m256i pix = _mm256_cvtepu8_epi32(_mm_set_epi64x(0,  *((uint64_t *)&lineIn[x + xmin]) ));
            __m256i mul = _mm256_mullo_epi32(mmk, pix);
            sss256 = _mm256_add_epi32(sss256, mul);
        }
        sss = _mm_add_epi32(
            _mm256_extractf128_si256(sss256, 0),
            _mm256_extractf128_si256(sss256, 1)
        );

#else

        sss = _mm_set1_epi32(1 << (PRECISION_BITS -1));
        for (; x < xmax - 3; x += 4) {
            __m128i pix, mmk, mul;
            __m128i source = _mm_loadu_si128((__m128i *) &lineIn[x + xmin]);
            __m128i ksource = _mm_loadu_si128((__m128i *) &k[x]);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            mmk = _mm_shuffle_epi32(ksource, _MM_SHUFFLE(0, 0, 0, 0));
            mul = _mm_mullo_epi32(pix, mmk);
            sss = _mm_add_epi32(sss, mul);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            mmk = _mm_shuffle_epi32(ksource, _MM_SHUFFLE(1, 1, 1, 1));
            mul = _mm_mullo_epi32(pix, mmk);
            sss = _mm_add_epi32(sss, mul);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            mmk = _mm_shuffle_epi32(ksource, _MM_SHUFFLE(2, 2, 2, 2));
            mul = _mm_mullo_epi32(pix, mmk);
            sss = _mm_add_epi32(sss, mul);

            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            mmk = _mm_shuffle_epi32(ksource, _MM_SHUFFLE(3, 3, 3, 3));
            mul = _mm_mullo_epi32(pix, mmk);
            sss = _mm_add_epi32(sss, mul);
        }

#endif

        for (; x < xmax; x ++) {
            __m128i pix = _mm_cvtepu8_epi32(*(__m128i *) &lineIn[x + xmin]);
            __m128i mmk = _mm_set1_epi32(k[x]);
            __m128i mul = _mm_mullo_epi32(pix, mmk);
            sss = _mm_add_epi32(sss, mul);
        }
        // sss = _mm_max_epi32(mmmin, _mm_min_epi32(mmmax, sss));
        sss = _mm_srai_epi32(sss, PRECISION_BITS);
        sss = _mm_packs_epi32(sss, sss);
        lineOut[xx] = _mm_cvtsi128_si32(_mm_packus_epi16(sss, sss));
    }
}

void
ImagingResampleVerticalConvolution8u(UINT32 *lineOut, Imaging imIn,
    int xmin, int xmax, int *k)
{
    int x;
    int xx = 0;
    int xsize = imIn->xsize;

#if defined(__AVX2__)

    for (; xx < xsize - 7; xx += 8) {
        __m256i pix, mmk;
        __m256i sss0 = _mm256_set1_epi32(1 << (PRECISION_BITS -1));
        __m256i sss1 = _mm256_set1_epi32(1 << (PRECISION_BITS -1));
        __m256i sss2 = _mm256_set1_epi32(1 << (PRECISION_BITS -1));
        __m256i sss3 = _mm256_set1_epi32(1 << (PRECISION_BITS -1));

        for (x = 0; x < xmax; x++) {
            __m256i source = _mm256_loadu_si256((__m256i *) &imIn->image32[x + xmin][xx]);
            mmk = _mm256_set1_epi32(k[x]);
            
            pix = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0,
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            sss0 = _mm256_add_epi32(sss0, _mm256_mullo_epi32(pix, mmk));
            
            pix = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4,
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            sss1 = _mm256_add_epi32(sss1, _mm256_mullo_epi32(pix, mmk));
            
            pix = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8,
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            sss2 = _mm256_add_epi32(sss2, _mm256_mullo_epi32(pix, mmk));
            
            pix = _mm256_shuffle_epi8(source, _mm256_set_epi8(
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12,
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            sss3 = _mm256_add_epi32(sss3, _mm256_mullo_epi32(pix, mmk));
        }
        sss0 = _mm256_srai_epi32(sss0, PRECISION_BITS);
        sss1 = _mm256_srai_epi32(sss1, PRECISION_BITS);
        sss2 = _mm256_srai_epi32(sss2, PRECISION_BITS);
        sss3 = _mm256_srai_epi32(sss3, PRECISION_BITS);

        sss0 = _mm256_packus_epi32(sss0, sss1);
        sss2 = _mm256_packus_epi32(sss2, sss3);
        sss0 = _mm256_packus_epi16(sss0, sss2);
        _mm256_storeu_si256((__m256i *) &lineOut[xx], sss0);
    }

#endif
    
    for (; xx < xsize - 3; xx += 4) {
        __m128i pix, mmk;
        __m128i sss0 = _mm_set1_epi32(1 << (PRECISION_BITS -1));
        __m128i sss1 = _mm_set1_epi32(1 << (PRECISION_BITS -1));
        __m128i sss2 = _mm_set1_epi32(1 << (PRECISION_BITS -1));
        __m128i sss3 = _mm_set1_epi32(1 << (PRECISION_BITS -1));

        for (x = 0; x < xmax; x++) {
            __m128i source = _mm_loadu_si128((__m128i *) &imIn->image32[x + xmin][xx]);
            mmk = _mm_set1_epi32(k[x]);
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,3, -1,-1,-1,2, -1,-1,-1,1, -1,-1,-1,0));
            sss0 = _mm_add_epi32(sss0, _mm_mullo_epi32(pix, mmk));
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,7, -1,-1,-1,6, -1,-1,-1,5, -1,-1,-1,4));
            sss1 = _mm_add_epi32(sss1, _mm_mullo_epi32(pix, mmk));
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,11, -1,-1,-1,10, -1,-1,-1,9, -1,-1,-1,8));
            sss2 = _mm_add_epi32(sss2, _mm_mullo_epi32(pix, mmk));
            
            pix = _mm_shuffle_epi8(source, _mm_set_epi8(
                -1,-1,-1,15, -1,-1,-1,14, -1,-1,-1,13, -1,-1,-1,12));
            sss3 = _mm_add_epi32(sss3, _mm_mullo_epi32(pix, mmk));
        }
        sss0 = _mm_srai_epi32(sss0, PRECISION_BITS);
        sss1 = _mm_srai_epi32(sss1, PRECISION_BITS);
        sss2 = _mm_srai_epi32(sss2, PRECISION_BITS);
        sss3 = _mm_srai_epi32(sss3, PRECISION_BITS);

        sss0 = _mm_packs_epi32(sss0, sss1);
        sss2 = _mm_packs_epi32(sss2, sss3);
        sss0 = _mm_packus_epi16(sss0, sss2);
        _mm_storeu_si128((__m128i *) &lineOut[xx], sss0);
    }

    for (; xx < imIn->xsize; xx++) {
        __m128i pix, mmk;
        __m128i sss = _mm_set1_epi32(1 << (PRECISION_BITS -1));
        for (x = 0; x < xmax; x++) {
            pix = _mm_cvtepu8_epi32(*(__m128i *) &imIn->image32[x + xmin][xx]);
            mmk = _mm_set1_epi32(k[x]);
            sss = _mm_add_epi32(sss, _mm_mullo_epi32(pix, mmk));
        }
        sss = _mm_srai_epi32(sss, PRECISION_BITS);
        sss = _mm_packs_epi32(sss, sss);
        lineOut[xx] = _mm_cvtsi128_si32(_mm_packus_epi16(sss, sss));
    }
}


Imaging
ImagingResampleHorizontal_8bpc(Imaging imIn, int xsize, struct filter *filterp)
{
    ImagingSectionCookie cookie;
    Imaging imOut;
    int ss0;
    int xx, yy, x, kmax, xmin, xmax;
    int *xbounds;
    int *k, *kk;
    double *prekk;


    kmax = ImagingPrecompute(imIn->xsize, xsize, filterp, &xbounds, &prekk);
    if ( ! kmax) {
        return (Imaging) ImagingError_MemoryError();
    }
    
    kk = malloc(xsize * kmax * sizeof(int));
    if ( ! kk) {
        free(xbounds);
        free(prekk);
        return (Imaging) ImagingError_MemoryError();
    }

    for (x = 0; x < xsize * kmax; x++) {
        kk[x] = (int) (0.5 + prekk[x] * (1 << PRECISION_BITS));
    }

    free(prekk);

    imOut = ImagingNew(imIn->mode, xsize, imIn->ysize);
    if ( ! imOut) {
        free(kk);
        free(xbounds);
        return NULL;
    }

    ImagingSectionEnter(&cookie);
    if (imIn->image8) {
        for (yy = 0; yy < imOut->ysize; yy++) {
            for (xx = 0; xx < xsize; xx++) {
                xmin = xbounds[xx * 2 + 0];
                xmax = xbounds[xx * 2 + 1];
                k = &kk[xx * kmax];
                ss0 = 1 << (PRECISION_BITS -1);
                for (x = 0; x < xmax; x++)
                    ss0 += ((UINT8) imIn->image8[yy][x + xmin]) * k[x];
                imOut->image8[yy][xx] = clip8(ss0);
            }
        }
    } else if (imIn->type == IMAGING_TYPE_UINT8) {
        for (yy = 0; yy < imOut->ysize; yy++) {
            ImagingResampleHorizontalConvolution8u(
                (UINT32 *) imOut->image32[yy],
                (UINT32 *) imIn->image32[yy],
                xsize, xbounds, kk, kmax
            );
        }
    }

    ImagingSectionLeave(&cookie);
    free(kk);
    free(xbounds);
    return imOut;
}


Imaging
ImagingResampleVertical_8bpc(Imaging imIn, int ysize, struct filter *filterp)
{
    ImagingSectionCookie cookie;
    Imaging imOut;
    int ss0;
    int xx, yy, y, kmax, ymin, ymax;
    int *xbounds;
    int *k, *kk;
    double *prekk;


    kmax = ImagingPrecompute(imIn->ysize, ysize, filterp, &xbounds, &prekk);
    if ( ! kmax) {
        return (Imaging) ImagingError_MemoryError();
    }
    
    kk = malloc(ysize * kmax * sizeof(int));
    if ( ! kk) {
        free(xbounds);
        free(prekk);
        return (Imaging) ImagingError_MemoryError();
    }

    for (y = 0; y < ysize * kmax; y++) {
        kk[y] = (int) (0.5 + prekk[y] * (1 << PRECISION_BITS));
    }

    free(prekk);

    imOut = ImagingNew(imIn->mode, imIn->xsize, ysize);
    if ( ! imOut) {
        free(kk);
        free(xbounds);
        return NULL;
    }

    ImagingSectionEnter(&cookie);
    if (imIn->image8) {
        for (yy = 0; yy < ysize; yy++) {
            k = &kk[yy * kmax];
            ymin = xbounds[yy * 2 + 0];
            ymax = xbounds[yy * 2 + 1];
            for (xx = 0; xx < imOut->xsize; xx++) {
                ss0 = 1 << (PRECISION_BITS -1);
                for (y = 0; y < ymax; y++)
                    ss0 += ((UINT8) imIn->image8[y + ymin][xx]) * k[y];
                imOut->image8[yy][xx] = clip8(ss0);
            }
        }
    } else if (imIn->type == IMAGING_TYPE_UINT8) {
        for (yy = 0; yy < ysize; yy++) {
            k = &kk[yy * kmax];
            ymin = xbounds[yy * 2 + 0];
            ymax = xbounds[yy * 2 + 1];
            ImagingResampleVerticalConvolution8u(
                (UINT32 *) imOut->image32[yy], imIn,
                ymin, ymax, k
            );
        }
    }

    ImagingSectionLeave(&cookie);
    free(kk);
    free(xbounds);
    return imOut;
}


Imaging
ImagingResampleHorizontal_32bpc(Imaging imIn, int xsize, struct filter *filterp)
{
    ImagingSectionCookie cookie;
    Imaging imOut;
    double ss;
    int xx, yy, x, kmax, xmin, xmax;
    int *xbounds;
    double *k, *kk;

    kmax = ImagingPrecompute(imIn->xsize, xsize, filterp, &xbounds, &kk);
    if ( ! kmax) {
        return (Imaging) ImagingError_MemoryError();
    }

    imOut = ImagingNew(imIn->mode, xsize, imIn->ysize);
    if ( ! imOut) {
        free(kk);
        free(xbounds);
        return NULL;
    }

    ImagingSectionEnter(&cookie);
    switch(imIn->type) {
        case IMAGING_TYPE_INT32:
            for (yy = 0; yy < imOut->ysize; yy++) {
                for (xx = 0; xx < xsize; xx++) {
                    xmin = xbounds[xx * 2 + 0];
                    xmax = xbounds[xx * 2 + 1];
                    k = &kk[xx * kmax];
                    ss = 0.0;
                    for (x = 0; x < xmax; x++)
                        ss += IMAGING_PIXEL_I(imIn, x + xmin, yy) * k[x];
                    IMAGING_PIXEL_I(imOut, xx, yy) = ROUND_UP(ss);
                }
            }
            break;

        case IMAGING_TYPE_FLOAT32:
            for (yy = 0; yy < imOut->ysize; yy++) {
                for (xx = 0; xx < xsize; xx++) {
                    xmin = xbounds[xx * 2 + 0];
                    xmax = xbounds[xx * 2 + 1];
                    k = &kk[xx * kmax];
                    ss = 0.0;
                    for (x = 0; x < xmax; x++)
                        ss += IMAGING_PIXEL_F(imIn, x + xmin, yy) * k[x];
                    IMAGING_PIXEL_F(imOut, xx, yy) = ss;
                }
            }
            break;
    }

    ImagingSectionLeave(&cookie);
    free(kk);
    free(xbounds);
    return imOut;
}


Imaging
ImagingResampleVertical_32bpc(Imaging imIn, int ysize, struct filter *filterp)
{
    ImagingSectionCookie cookie;
    Imaging imOut;
    double ss;
    int xx, yy, y, kmax, ymin, ymax;
    int *xbounds;
    double *k, *kk;

    kmax = ImagingPrecompute(imIn->ysize, ysize, filterp, &xbounds, &kk);
    if ( ! kmax) {
        return (Imaging) ImagingError_MemoryError();
    }

    imOut = ImagingNew(imIn->mode, imIn->xsize, ysize);
    if ( ! imOut) {
        free(kk);
        free(xbounds);
        return NULL;
    }

    ImagingSectionEnter(&cookie);
    switch(imIn->type) {
        case IMAGING_TYPE_INT32:
            for (yy = 0; yy < ysize; yy++) {
                ymin = xbounds[yy * 2 + 0];
                ymax = xbounds[yy * 2 + 1];
                k = &kk[yy * kmax];
                for (xx = 0; xx < imOut->xsize; xx++) {
                    ss = 0.0;
                    for (y = 0; y < ymax; y++)
                        ss += IMAGING_PIXEL_I(imIn, xx, y + ymin) * k[y];
                    IMAGING_PIXEL_I(imOut, xx, yy) = ROUND_UP(ss);
                }
            }
            break;

        case IMAGING_TYPE_FLOAT32:
            for (yy = 0; yy < ysize; yy++) {
                ymin = xbounds[yy * 2 + 0];
                ymax = xbounds[yy * 2 + 1];
                k = &kk[yy * kmax];
                for (xx = 0; xx < imOut->xsize; xx++) {
                    ss = 0.0;
                    for (y = 0; y < ymax; y++)
                        ss += IMAGING_PIXEL_F(imIn, xx, y + ymin) * k[y];
                    IMAGING_PIXEL_F(imOut, xx, yy) = ss;
                }
            }
            break;
    }

    ImagingSectionLeave(&cookie);
    free(kk);
    free(xbounds);
    return imOut;
}


Imaging
ImagingResample(Imaging imIn, int xsize, int ysize, int filter)
{
    Imaging imTemp = NULL;
    Imaging imOut = NULL;
    struct filter *filterp;
    Imaging (*ResampleHorizontal)(Imaging imIn, int xsize, struct filter *filterp);
    Imaging (*ResampleVertical)(Imaging imIn, int xsize, struct filter *filterp);

    if (strcmp(imIn->mode, "P") == 0 || strcmp(imIn->mode, "1") == 0)
        return (Imaging) ImagingError_ModeError();

    if (imIn->type == IMAGING_TYPE_SPECIAL) {
        return (Imaging) ImagingError_ModeError();
    } else if (imIn->image8) {
        ResampleHorizontal = ImagingResampleHorizontal_8bpc;
        ResampleVertical = ImagingResampleVertical_8bpc;
    } else {
        switch(imIn->type) {
            case IMAGING_TYPE_UINT8:
                ResampleHorizontal = ImagingResampleHorizontal_8bpc;
                ResampleVertical = ImagingResampleVertical_8bpc;
                break;
            case IMAGING_TYPE_INT32:
            case IMAGING_TYPE_FLOAT32:
                ResampleHorizontal = ImagingResampleHorizontal_32bpc;
                ResampleVertical = ImagingResampleVertical_32bpc;
                break;
            default:
                return (Imaging) ImagingError_ModeError();
        }
    }

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

    /* two-pass resize, first pass */
    if (imIn->xsize != xsize) {
        imTemp = ResampleHorizontal(imIn, xsize, filterp);
        if ( ! imTemp)
            return NULL;
        imOut = imIn = imTemp;
    }

    /* second pass */
    if (imIn->ysize != ysize) {
        /* imIn can be the original image or horizontally resampled one */
        imOut = ResampleVertical(imIn, ysize, filterp);
        /* it's safe to call ImagingDelete with empty value
           if there was no previous step. */
        ImagingDelete(imTemp);
        if ( ! imOut)
            return NULL;
    }

    /* none of the previous steps are performed, copying */
    if ( ! imOut) {
        imOut = ImagingCopy(imIn);
    }

    return imOut;
}
