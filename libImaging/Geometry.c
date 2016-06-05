#include "Imaging.h"

/* For large images rotation is an inefficient operation in terms of CPU cache.
   One row in the source image affects each column in destination.
   Rotating in chunks that fit in the cache can speed up rotation
   8x on a modern CPU. A chunk size of 128 requires only 65k and is large enough
   that the overhead from the extra loops are not apparent. */
#define ROTATE_CHUNK 128

#define COORD(v) ((v) < 0.0 ? -1 : ((int)(v)))
#define FLOOR(v) ((v) < 0.0 ? ((int)floor(v)) : ((int)(v)))

/* -------------------------------------------------------------------- */
/* Transpose operations                                                 */

Imaging
ImagingFlipLeftRight(Imaging imOut, Imaging imIn)
{
    ImagingSectionCookie cookie;
    int x, y, xr;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();
    if (imIn->xsize != imOut->xsize || imIn->ysize != imOut->ysize)
        return (Imaging) ImagingError_Mismatch();

    ImagingCopyInfo(imOut, imIn);

#define FLIP_HORIZ(image)\
    for (y = 0; y < imIn->ysize; y++) {\
        xr = imIn->xsize-1;\
        for (x = 0; x < imIn->xsize; x++, xr--)\
            imOut->image[y][x] = imIn->image[y][xr];\
    }

    ImagingSectionEnter(&cookie);

    if (imIn->image8)
        FLIP_HORIZ(image8)
    else
        FLIP_HORIZ(image32)

    ImagingSectionLeave(&cookie);

#undef FLIP_HORIZ

    return imOut;
}


Imaging
ImagingFlipTopBottom(Imaging imOut, Imaging imIn)
{
    ImagingSectionCookie cookie;
    int y, yr;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();
    if (imIn->xsize != imOut->xsize || imIn->ysize != imOut->ysize)
        return (Imaging) ImagingError_Mismatch();

    ImagingCopyInfo(imOut, imIn);

    ImagingSectionEnter(&cookie);

    yr = imIn->ysize - 1;
    for (y = 0; y < imIn->ysize; y++, yr--)
        memcpy(imOut->image[yr], imIn->image[y], imIn->linesize);

    ImagingSectionLeave(&cookie);

    return imOut;
}


Imaging
ImagingRotate90(Imaging imOut, Imaging imIn)
{
    ImagingSectionCookie cookie;
    int x, y, xx, yy, xr, xxsize, yysize;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();
    if (imIn->xsize != imOut->ysize || imIn->ysize != imOut->xsize)
        return (Imaging) ImagingError_Mismatch();

    ImagingCopyInfo(imOut, imIn);

#define ROTATE_90(image) \
    for (y = 0; y < imIn->ysize; y += ROTATE_CHUNK) { \
        for (x = 0; x < imIn->xsize; x += ROTATE_CHUNK) { \
            yysize = y + ROTATE_CHUNK < imIn->ysize ? y + ROTATE_CHUNK : imIn->ysize; \
            xxsize = x + ROTATE_CHUNK < imIn->xsize ? x + ROTATE_CHUNK : imIn->xsize; \
            for (yy = y; yy < yysize; yy++) { \
                xr = imIn->xsize - 1 - x; \
                for (xx = x; xx < xxsize; xx++, xr--) { \
                    imOut->image[xr][yy] = imIn->image[yy][xx]; \
                } \
            } \
        } \
    }

    ImagingSectionEnter(&cookie);

    if (imIn->image8)
        ROTATE_90(image8)
    else
        ROTATE_90(image32)

    ImagingSectionLeave(&cookie);

#undef ROTATE_90

    return imOut;
}


Imaging
ImagingTranspose(Imaging imOut, Imaging imIn)
{
    ImagingSectionCookie cookie;
    int x, y, xx, yy, xxsize, yysize;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();
    if (imIn->xsize != imOut->ysize || imIn->ysize != imOut->xsize)
        return (Imaging) ImagingError_Mismatch();

    ImagingCopyInfo(imOut, imIn);

#define TRANSPOSE(image) \
    for (y = 0; y < imIn->ysize; y += ROTATE_CHUNK) { \
        for (x = 0; x < imIn->xsize; x += ROTATE_CHUNK) { \
            yysize = y + ROTATE_CHUNK < imIn->ysize ? y + ROTATE_CHUNK : imIn->ysize; \
            xxsize = x + ROTATE_CHUNK < imIn->xsize ? x + ROTATE_CHUNK : imIn->xsize; \
            for (yy = y; yy < yysize; yy++) { \
                for (xx = x; xx < xxsize; xx++) { \
                    imOut->image[xx][yy] = imIn->image[yy][xx]; \
                } \
            } \
        } \
    }

    ImagingSectionEnter(&cookie);
    
    if (imIn->image8)
        TRANSPOSE(image8)
    else
        TRANSPOSE(image32)

    ImagingSectionLeave(&cookie);

#undef TRANSPOSE

    return imOut;
}


Imaging
ImagingRotate180(Imaging imOut, Imaging imIn)
{
    ImagingSectionCookie cookie;
    int x, y, xr, yr;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();
    if (imIn->xsize != imOut->xsize || imIn->ysize != imOut->ysize)
        return (Imaging) ImagingError_Mismatch();

    ImagingCopyInfo(imOut, imIn);

#define ROTATE_180(image)\
    for (y = 0; y < imIn->ysize; y++, yr--) {\
        xr = imIn->xsize-1;\
        for (x = 0; x < imIn->xsize; x++, xr--)\
            imOut->image[y][x] = imIn->image[yr][xr];\
    }

    ImagingSectionEnter(&cookie);

    yr = imIn->ysize-1;
    if (imIn->image8)
        ROTATE_180(image8)
    else
        ROTATE_180(image32)

    ImagingSectionLeave(&cookie);

#undef ROTATE_180

    return imOut;
}


Imaging
ImagingRotate270(Imaging imOut, Imaging imIn)
{
    ImagingSectionCookie cookie;
    int x, y, xx, yy, yr, xxsize, yysize;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();
    if (imIn->xsize != imOut->ysize || imIn->ysize != imOut->xsize)
        return (Imaging) ImagingError_Mismatch();

    ImagingCopyInfo(imOut, imIn);

#define ROTATE_270(image) \
    for (y = 0; y < imIn->ysize; y += ROTATE_CHUNK) { \
        for (x = 0; x < imIn->xsize; x += ROTATE_CHUNK) { \
            yysize = y + ROTATE_CHUNK < imIn->ysize ? y + ROTATE_CHUNK : imIn->ysize; \
            xxsize = x + ROTATE_CHUNK < imIn->xsize ? x + ROTATE_CHUNK : imIn->xsize; \
            yr = imIn->ysize - 1 - y; \
            for (yy = y; yy < yysize; yy++, yr--) { \
                for (xx = x; xx < xxsize; xx++) { \
                    imOut->image[xx][yr] = imIn->image[yy][xx]; \
                } \
            } \
        } \
    }

    ImagingSectionEnter(&cookie);

    if (imIn->image8)
        ROTATE_270(image8)
    else
        ROTATE_270(image32)

    ImagingSectionLeave(&cookie);

#undef ROTATE_270

    return imOut;
}


/* -------------------------------------------------------------------- */
/* Transforms                                                           */

/* transform primitives (ImagingTransformMap) */

static int
affine_transform(double* xin, double* yin, int x, int y, void* data)
{
    /* full moon tonight.  your compiler will generate bogus code
       for simple expressions, unless you reorganize the code, or
       install Service Pack 3 */

    double* a = (double*) data;
    double a0 = a[0]; double a1 = a[1]; double a2 = a[2];
    double a3 = a[3]; double a4 = a[4]; double a5 = a[5];

    xin[0] = a0*x + a1*y + a2;
    yin[0] = a3*x + a4*y + a5;

    return 1;
}

static int
perspective_transform(double* xin, double* yin, int x, int y, void* data)
{
    double* a = (double*) data;
    double a0 = a[0]; double a1 = a[1]; double a2 = a[2];
    double a3 = a[3]; double a4 = a[4]; double a5 = a[5];
    double a6 = a[6]; double a7 = a[7];

    xin[0] = (a0*x + a1*y + a2) / (a6*x + a7*y + 1);
    yin[0] = (a3*x + a4*y + a5) / (a6*x + a7*y + 1);

    return 1;
}

static int
quad_transform(double* xin, double* yin, int x, int y, void* data)
{
    /* quad warp: map quadrilateral to rectangle */

    double* a = (double*) data;
    double a0 = a[0]; double a1 = a[1]; double a2 = a[2]; double a3 = a[3];
    double a4 = a[4]; double a5 = a[5]; double a6 = a[6]; double a7 = a[7];

    xin[0] = a0 + a1*x + a2*y + a3*x*y;
    yin[0] = a4 + a5*x + a6*y + a7*x*y;

    return 1;
}

/* transform filters (ImagingTransformFilter) */

static int
nearest_filter8(void* out, Imaging im, double xin, double yin)
{
    int x = COORD(xin);
    int y = COORD(yin);
    if (x < 0 || x >= im->xsize || y < 0 || y >= im->ysize)
        return 0;
    ((UINT8*)out)[0] = im->image8[y][x];
    return 1;
}

static int
nearest_filter16(void* out, Imaging im, double xin, double yin)
{
    int x = COORD(xin);
    int y = COORD(yin);
    if (x < 0 || x >= im->xsize || y < 0 || y >= im->ysize)
        return 0;
    ((INT16*)out)[0] = ((INT16*)(im->image8[y]))[x];
    return 1;
}

static int
nearest_filter32(void* out, Imaging im, double xin, double yin)
{
    int x = COORD(xin);
    int y = COORD(yin);
    if (x < 0 || x >= im->xsize || y < 0 || y >= im->ysize)
        return 0;
    ((INT32*)out)[0] = im->image32[y][x];
    return 1;
}

#define XCLIP(im, x) ( ((x) < 0) ? 0 : ((x) < im->xsize) ? (x) : im->xsize-1 )
#define YCLIP(im, y) ( ((y) < 0) ? 0 : ((y) < im->ysize) ? (y) : im->ysize-1 )

#define BILINEAR(v, a, b, d)\
    (v = (a) + ( (b) - (a) ) * (d))

#define BILINEAR_HEAD(type)\
    int x, y;\
    int x0, x1;\
    double v1, v2;\
    double dx, dy;\
    type* in;\
    if (xin < 0.0 || xin >= im->xsize || yin < 0.0 || yin >= im->ysize)\
        return 0;\
    xin -= 0.5;\
    yin -= 0.5;\
    x = FLOOR(xin);\
    y = FLOOR(yin);\
    dx = xin - x;\
    dy = yin - y;

#define BILINEAR_BODY(type, image, step, offset) {\
    in = (type*) ((image)[YCLIP(im, y)] + offset);\
    x0 = XCLIP(im, x+0)*step;\
    x1 = XCLIP(im, x+1)*step;\
    BILINEAR(v1, in[x0], in[x1], dx);\
    if (y+1 >= 0 && y+1 < im->ysize) {\
        in = (type*) ((image)[y+1] + offset);\
        BILINEAR(v2, in[x0], in[x1], dx);\
    } else\
        v2 = v1;\
    BILINEAR(v1, v1, v2, dy);\
}

static int
bilinear_filter8(void* out, Imaging im, double xin, double yin)
{
    BILINEAR_HEAD(UINT8);
    BILINEAR_BODY(UINT8, im->image8, 1, 0);
    ((UINT8*)out)[0] = (UINT8) v1;
    return 1;
}

static int
bilinear_filter32I(void* out, Imaging im, double xin, double yin)
{
    BILINEAR_HEAD(INT32);
    BILINEAR_BODY(INT32, im->image32, 1, 0);
    ((INT32*)out)[0] = (INT32) v1;
    return 1;
}

static int
bilinear_filter32F(void* out, Imaging im, double xin, double yin)
{
    BILINEAR_HEAD(FLOAT32);
    BILINEAR_BODY(FLOAT32, im->image32, 1, 0);
    ((FLOAT32*)out)[0] = (FLOAT32) v1;
    return 1;
}

static int
bilinear_filter32LA(void* out, Imaging im, double xin, double yin)
{
    BILINEAR_HEAD(UINT8);
    BILINEAR_BODY(UINT8, im->image, 4, 0);
    ((UINT8*)out)[0] = (UINT8) v1;
    ((UINT8*)out)[1] = (UINT8) v1;
    ((UINT8*)out)[2] = (UINT8) v1;
    BILINEAR_BODY(UINT8, im->image, 4, 3);
    ((UINT8*)out)[3] = (UINT8) v1;
    return 1;
}

static int
bilinear_filter32RGB(void* out, Imaging im, double xin, double yin)
{
    int b;
    BILINEAR_HEAD(UINT8);
    for (b = 0; b < im->bands; b++) {
        BILINEAR_BODY(UINT8, im->image, 4, b);
        ((UINT8*)out)[b] = (UINT8) v1;
    }
    return 1;
}

#undef BILINEAR
#undef BILINEAR_HEAD
#undef BILINEAR_BODY

#define BICUBIC(v, v1, v2, v3, v4, d) {\
    double p1 = v2;\
    double p2 = -v1 + v3;\
    double p3 = 2*(v1 - v2) + v3 - v4;\
    double p4 = -v1 + v2 - v3 + v4;\
    v = p1 + (d)*(p2 + (d)*(p3 + (d)*p4));\
}

#define BICUBIC_HEAD(type)\
    int x = FLOOR(xin);\
    int y = FLOOR(yin);\
    int x0, x1, x2, x3;\
    double v1, v2, v3, v4;\
    double dx, dy;\
    type* in;\
    if (xin < 0.0 || xin >= im->xsize || yin < 0.0 || yin >= im->ysize)\
        return 0;\
    xin -= 0.5;\
    yin -= 0.5;\
    x = FLOOR(xin);\
    y = FLOOR(yin);\
    dx = xin - x;\
    dy = yin - y;\
    x--; y--;

#define BICUBIC_BODY(type, image, step, offset) {\
    in = (type*) ((image)[YCLIP(im, y)] + offset);\
    x0 = XCLIP(im, x+0)*step;\
    x1 = XCLIP(im, x+1)*step;\
    x2 = XCLIP(im, x+2)*step;\
    x3 = XCLIP(im, x+3)*step;\
    BICUBIC(v1, in[x0], in[x1], in[x2], in[x3], dx);\
    if (y+1 >= 0 && y+1 < im->ysize) {\
        in = (type*) ((image)[y+1] + offset);\
        BICUBIC(v2, in[x0], in[x1], in[x2], in[x3], dx);\
    } else\
        v2 = v1;\
    if (y+2 >= 0 && y+2 < im->ysize) {\
        in = (type*) ((image)[y+2] + offset);\
        BICUBIC(v3, in[x0], in[x1], in[x2], in[x3], dx);\
    } else\
        v3 = v2;\
    if (y+3 >= 0 && y+3 < im->ysize) {\
        in = (type*) ((image)[y+3] + offset);\
        BICUBIC(v4, in[x0], in[x1], in[x2], in[x3], dx);\
    } else\
        v4 = v3;\
    BICUBIC(v1, v1, v2, v3, v4, dy);\
}


static int
bicubic_filter8(void* out, Imaging im, double xin, double yin)
{
    BICUBIC_HEAD(UINT8);
    BICUBIC_BODY(UINT8, im->image8, 1, 0);
    if (v1 <= 0.0)
        ((UINT8*)out)[0] = 0;
    else if (v1 >= 255.0)
        ((UINT8*)out)[0] = 255;
    else
        ((UINT8*)out)[0] = (UINT8) v1;
    return 1;
}

static int
bicubic_filter32I(void* out, Imaging im, double xin, double yin)
{
    BICUBIC_HEAD(INT32);
    BICUBIC_BODY(INT32, im->image32, 1, 0);
    ((INT32*)out)[0] = (INT32) v1;
    return 1;
}

static int
bicubic_filter32F(void* out, Imaging im, double xin, double yin)
{
    BICUBIC_HEAD(FLOAT32);
    BICUBIC_BODY(FLOAT32, im->image32, 1, 0);
    ((FLOAT32*)out)[0] = (FLOAT32) v1;
    return 1;
}

static int
bicubic_filter32LA(void* out, Imaging im, double xin, double yin)
{
    BICUBIC_HEAD(UINT8);
    BICUBIC_BODY(UINT8, im->image, 4, 0);
    if (v1 <= 0.0) {
        ((UINT8*)out)[0] = 0;
        ((UINT8*)out)[1] = 0;
        ((UINT8*)out)[2] = 0;
    } else if (v1 >= 255.0) {
        ((UINT8*)out)[0] = 255;
        ((UINT8*)out)[1] = 255;
        ((UINT8*)out)[2] = 255;
    } else {
        ((UINT8*)out)[0] = (UINT8) v1;
        ((UINT8*)out)[1] = (UINT8) v1;
        ((UINT8*)out)[2] = (UINT8) v1;
    }
    BICUBIC_BODY(UINT8, im->image, 4, 3);
    if (v1 <= 0.0)
        ((UINT8*)out)[3] = 0;
    else if (v1 >= 255.0)
        ((UINT8*)out)[3] = 255;
    else
        ((UINT8*)out)[3] = (UINT8) v1;
    return 1;
}

static int
bicubic_filter32RGB(void* out, Imaging im, double xin, double yin)
{
    int b;
    BICUBIC_HEAD(UINT8);
    for (b = 0; b < im->bands; b++) {
        BICUBIC_BODY(UINT8, im->image, 4, b);
        if (v1 <= 0.0)
            ((UINT8*)out)[b] = 0;
        else if (v1 >= 255.0)
            ((UINT8*)out)[b] = 255;
        else
            ((UINT8*)out)[b] = (UINT8) v1;
    }
    return 1;
}

#undef BICUBIC
#undef BICUBIC_HEAD
#undef BICUBIC_BODY

static ImagingTransformFilter
getfilter(Imaging im, int filterid)
{
    switch (filterid) {
    case IMAGING_TRANSFORM_NEAREST:
        if (im->image8)
            switch (im->type) {
            case IMAGING_TYPE_UINT8:
                return nearest_filter8;
            case IMAGING_TYPE_SPECIAL:
                switch (im->pixelsize) {
                case 1:
                    return nearest_filter8;
                case 2:
                    return nearest_filter16;
                case 4:
                    return nearest_filter32;
                }
            }
        else
            return nearest_filter32;
        break;
    case IMAGING_TRANSFORM_BILINEAR:
        if (im->image8)
            return bilinear_filter8;
        else if (im->image32) {
            switch (im->type) {
            case IMAGING_TYPE_UINT8:
                if (im->bands == 2)
                    return bilinear_filter32LA;
                else
                    return bilinear_filter32RGB;
            case IMAGING_TYPE_INT32:
                return bilinear_filter32I;
            case IMAGING_TYPE_FLOAT32:
                return bilinear_filter32F;
            }
        }
        break;
    case IMAGING_TRANSFORM_BICUBIC:
        if (im->image8)
            return bicubic_filter8;
        else if (im->image32) {
            switch (im->type) {
            case IMAGING_TYPE_UINT8:
                if (im->bands == 2)
                    return bicubic_filter32LA;
                else
                    return bicubic_filter32RGB;
            case IMAGING_TYPE_INT32:
                return bicubic_filter32I;
            case IMAGING_TYPE_FLOAT32:
                return bicubic_filter32F;
            }
        }
        break;
    }
    /* no such filter */
    return NULL;
}

/* transformation engines */

Imaging
ImagingGenericTransform(
    Imaging imOut, Imaging imIn, int x0, int y0, int x1, int y1,
    ImagingTransformMap transform, void* transform_data,
    int filterid, int fill)
{
    /* slow generic transformation.  use ImagingTransformAffine or
       ImagingScaleAffine where possible. */

    ImagingSectionCookie cookie;
    int x, y;
    char *out;
    double xx, yy;

    ImagingTransformFilter filter = getfilter(imIn, filterid);
    if (!filter)
        return (Imaging) ImagingError_ValueError("bad filter number");

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();

    ImagingCopyInfo(imOut, imIn);

    ImagingSectionEnter(&cookie);

    if (x0 < 0)
        x0 = 0;
    if (y0 < 0)
        y0 = 0;
    if (x1 > imOut->xsize)
        x1 = imOut->xsize;
    if (y1 > imOut->ysize)
        y1 = imOut->ysize;

    for (y = y0; y < y1; y++) {
        out = imOut->image[y] + x0*imOut->pixelsize;
        for (x = x0; x < x1; x++) {
            if ( ! transform(&xx, &yy, x-x0, y-y0, transform_data) ||
                 ! filter(out, imIn, xx, yy)) {
                if (fill)
                    memset(out, 0, imOut->pixelsize);
            }
            out += imOut->pixelsize;
        }
    }

    ImagingSectionLeave(&cookie);

    return imOut;
}

static Imaging
ImagingScaleAffine(Imaging imOut, Imaging imIn,
                   int x0, int y0, int x1, int y1,
                   double a[6], int fill)
{
    /* scale, nearest neighbour resampling */

    ImagingSectionCookie cookie;
    int x, y;
    int xin;
    double xo, yo;
    int xmin, xmax;
    int *xintab;

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();

    ImagingCopyInfo(imOut, imIn);

    if (x0 < 0)
        x0 = 0;
    if (y0 < 0)
        y0 = 0;
    if (x1 > imOut->xsize)
        x1 = imOut->xsize;
    if (y1 > imOut->ysize)
        y1 = imOut->ysize;

    xintab = (int*) malloc(imOut->xsize * sizeof(int));
    if (!xintab) {
        ImagingDelete(imOut);
        return (Imaging) ImagingError_MemoryError();
    }

    xo = a[2];
    yo = a[5];

    xmin = x1;
    xmax = x0;

    /* Pretabulate horizontal pixel positions */
    for (x = x0; x < x1; x++) {
        xin = COORD(xo);
        if (xin >= 0 && xin < (int) imIn->xsize) {
            xmax = x+1;
            if (x < xmin)
                xmin = x;
            xintab[x] = xin;
        }
        xo += a[0];
    }

#define AFFINE_SCALE(pixel, image)\
    for (y = y0; y < y1; y++) {\
        int yi = COORD(yo);\
        pixel *in, *out;\
        out = imOut->image[y];\
        if (fill && x1 > x0)\
            memset(out+x0, 0, (x1-x0)*sizeof(pixel));\
        if (yi >= 0 && yi < imIn->ysize) {\
            in = imIn->image[yi];\
            for (x = xmin; x < xmax; x++)\
                out[x] = in[xintab[x]];\
        }\
        yo += a[4];\
    }

    ImagingSectionEnter(&cookie);

    if (imIn->image8) {
        AFFINE_SCALE(UINT8, image8);
    } else {
        AFFINE_SCALE(INT32, image32);
    }

    ImagingSectionLeave(&cookie);

#undef AFFINE_SCALE

    free(xintab);

    return imOut;
}

static inline int
check_fixed(double a[6], int x, int y)
{
    return (fabs(x*a[0] + y*a[1] + a[2]) < 32768.0 &&
            fabs(x*a[3] + y*a[4] + a[5]) < 32768.0);
}

static inline Imaging
affine_fixed(Imaging imOut, Imaging imIn,
             int x0, int y0, int x1, int y1,
             double a[6], int filterid, int fill)
{
    /* affine transform, nearest neighbour resampling, fixed point
       arithmetics */

    ImagingSectionCookie cookie;
    int x, y;
    int xin, yin;
    int xsize, ysize;
    int xx, yy;
    int a0, a1, a2, a3, a4, a5;

    ImagingCopyInfo(imOut, imIn);

    xsize = (int) imIn->xsize;
    ysize = (int) imIn->ysize;

/* use 16.16 fixed point arithmetics */
#define FIX(v) FLOOR((v)*65536.0 + 0.5)

    a0 = FIX(a[0]); a1 = FIX(a[1]); a2 = FIX(a[2]);
    a3 = FIX(a[3]); a4 = FIX(a[4]); a5 = FIX(a[5]);

#undef FIX

#define AFFINE_TRANSFORM_FIXED(pixel, image)\
    for (y = y0; y < y1; y++) {\
        pixel *out;\
        xx = a2;\
        yy = a5;\
        out = imOut->image[y];\
        if (fill && x1 > x0)\
            memset(out+x0, 0, (x1-x0)*sizeof(pixel));\
        for (x = x0; x < x1; x++, out++) {\
            xin = xx >> 16;\
            if (xin >= 0 && xin < xsize) {\
                yin = yy >> 16;\
                if (yin >= 0 && yin < ysize)\
                    *out = imIn->image[yin][xin];\
            }\
            xx += a0;\
            yy += a3;\
        }\
        a2 += a1;\
        a5 += a4;\
    }

    ImagingSectionEnter(&cookie);

    if (imIn->image8)
        AFFINE_TRANSFORM_FIXED(UINT8, image8)
    else
        AFFINE_TRANSFORM_FIXED(INT32, image32)

    ImagingSectionLeave(&cookie);

#undef AFFINE_TRANSFORM_FIXED

    return imOut;
}

Imaging
ImagingTransformAffine(Imaging imOut, Imaging imIn,
                       int x0, int y0, int x1, int y1,
                       double a[6], int filterid, int fill)
{
    /* affine transform, nearest neighbour resampling, floating point
       arithmetics*/

    ImagingSectionCookie cookie;
    int x, y;
    int xin, yin;
    int xsize, ysize;
    double xx, yy;
    double xo, yo;

    if (filterid || imIn->type == IMAGING_TYPE_SPECIAL) {
        return ImagingGenericTransform(
            imOut, imIn,
            x0, y0, x1, y1,
            affine_transform, a,
            filterid, fill);
    }

    if (a[1] == 0 && a[3] == 0)
        /* Scaling */
        return ImagingScaleAffine(imOut, imIn, x0, y0, x1, y1, a, fill);

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();

    if (x0 < 0)
        x0 = 0;
    if (y0 < 0)
        y0 = 0;
    if (x1 > imOut->xsize)
        x1 = imOut->xsize;
    if (y1 > imOut->ysize)
        y1 = imOut->ysize;

    /* translate all four corners to check if they are within the
       range that can be represented by the fixed point arithmetics */

    if (check_fixed(a, 0, 0) && check_fixed(a, x1-x0, y1-y0) &&
        check_fixed(a, 0, y1-y0) && check_fixed(a, x1-x0, 0))
        return affine_fixed(imOut, imIn, x0, y0, x1, y1, a, filterid, fill);

    /* FIXME: cannot really think of any reasonable case when the
       following code is used.  maybe we should fall back on the slow
       generic transform engine in this case? */
    
    ImagingCopyInfo(imOut, imIn);

    xsize = (int) imIn->xsize;
    ysize = (int) imIn->ysize;

    xo = a[2];
    yo = a[5];

#define AFFINE_TRANSFORM(pixel, image)\
    for (y = y0; y < y1; y++) {\
        pixel *out;\
        xx = xo;\
        yy = yo;\
        out = imOut->image[y];\
        if (fill && x1 > x0)\
            memset(out+x0, 0, (x1-x0)*sizeof(pixel));\
        for (x = x0; x < x1; x++, out++) {\
            xin = COORD(xx);\
            if (xin >= 0 && xin < xsize) {\
                yin = COORD(yy);\
                if (yin >= 0 && yin < ysize)\
                    *out = imIn->image[yin][xin];\
            }\
            xx += a[0];\
            yy += a[3];\
        }\
        xo += a[1];\
        yo += a[4];\
    }

    ImagingSectionEnter(&cookie);

    if (imIn->image8)
        AFFINE_TRANSFORM(UINT8, image8)
    else
        AFFINE_TRANSFORM(INT32, image32)

    ImagingSectionLeave(&cookie);

#undef AFFINE_TRANSFORM

    return imOut;
}

typedef struct point {
    float x;
    float y;
} point;


point inline
affine_transform2(int x, int y, double a[6]) {
    point r;
    r.x = a[0] * x + a[1] * y + a[2];
    r.y = a[3] * x + a[4] * y + a[5];
    return r;
}


struct filter {
    double (*filter)(double x);
    int support;
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

static struct filter LANCZOS = { lanczos_filter, 3 };
static struct filter BICUBIC = { bicubic_filter, 2 };
static struct filter BILINEAR = { bilinear_filter, 1 };


void
invert_affine_matrix(double m[6], double inv[6])
{
    double a = m[0], b = m[1], c = m[2];
    double d = m[3], e = m[4], f = m[5];
    double det = a*e - b*d;
    inv[0] =  e / det;
    inv[1] = -b / det;
    inv[3] = -d / det;
    inv[4] =  a / det;
    inv[2] = (b*f - e*c) / det;
    inv[5] = (c*d - a*f) / det;
}

UINT8 clip8(float inf)
{
    int in = (int) inf;
    if (in >= 255)
        return 255;
    if (in <= 0)
        return 0;
    return (UINT8) in;
}


Imaging
ImagingTransformAffineSlow(Imaging imOut, Imaging imIn,
                           int x0, int y0, int x1, int y1,
                           double a[6], int filterid, int fill)
{
    ImagingSectionCookie cookie;
    int x, y, xx, yy;
    float ss0, ss1, ss2, ss3;
    float xfrom, xto, yfrom, yto;
    float w, l;
    point pcenter, ptop, pbottom, pleft, pright;
    double inva[6];

    ImagingTransformFilter filter = getfilter(imIn, filterid);
    if (!filter)
        return (Imaging) ImagingError_ValueError("bad filter number");

    if (!imOut || !imIn || strcmp(imIn->mode, imOut->mode) != 0)
        return (Imaging) ImagingError_ModeError();

    ImagingCopyInfo(imOut, imIn);

    ImagingSectionEnter(&cookie);

    invert_affine_matrix(a, inva);

    if (x0 < 0)
        x0 = 0;
    if (y0 < 0)
        y0 = 0;
    if (x1 > imOut->xsize)
        x1 = imOut->xsize;
    if (y1 > imOut->ysize)
        y1 = imOut->ysize;

    for (y = y0; y < y1; y++) {
        for (x = x0; x < x1; x++) {
            ptop = affine_transform2(x, y + LANCZOS.support, a);
            pright = affine_transform2(x + LANCZOS.support, y, a);
            pbottom = affine_transform2(x, y - LANCZOS.support, a);
            pleft = affine_transform2(x - LANCZOS.support, y, a);

            xfrom = xto = ptop.x;
            if (xfrom > pright.x) xfrom = pright.x;
            if (xto < pright.x) xto = pright.x;
            if (xfrom > pbottom.x) xfrom = pbottom.x;
            if (xto < pbottom.x) xto = pbottom.x;
            if (xfrom > pleft.x) xfrom = pleft.x;
            if (xto < pleft.x) xto = pleft.x;
            yfrom = yto = ptop.y;
            if (yfrom > pright.y) yfrom = pright.y;
            if (yto < pright.y) yto = pright.y;
            if (yfrom > pbottom.y) yfrom = pbottom.y;
            if (yto < pbottom.y) yto = pbottom.y;
            if (yfrom > pleft.y) yfrom = pleft.y;
            if (yto < pleft.y) yto = pleft.y;

            if (yfrom < 0) yfrom = 0;
            if (yto > imIn->ysize) yto = imIn->ysize;
            if (xfrom < 0) xfrom = 0;
            if (xto > imIn->xsize) xto = imIn->xsize;

            float ww = 0;
            ss0 = ss1 = ss2 = ss3 = 0;
            for (yy = (int) yfrom; yy < yto; yy++) {
                for (xx = (int) xfrom; xx < xto; xx++) {
                    pcenter = affine_transform2(xx, yy, inva);
                    l = sqrt((pcenter.x - x) * (pcenter.x - x) +
                        (pcenter.y - y) * (pcenter.y - y));
                    if (l > LANCZOS.support)
                        continue;
                    w = LANCZOS.filter(l);

                    ss0 += ((UINT8) imIn->image[yy][xx*4 + 0]) * w;
                    ss1 += ((UINT8) imIn->image[yy][xx*4 + 1]) * w;
                    ss2 += ((UINT8) imIn->image[yy][xx*4 + 2]) * w;
                    ss3 += ((UINT8) imIn->image[yy][xx*4 + 3]) * w;
                    ww += w;
                }
            }
            if (ww) {
                imOut->image[y][x*4 + 0] = clip8(ss0 / ww);
                imOut->image[y][x*4 + 1] = clip8(ss1 / ww);
                imOut->image[y][x*4 + 2] = clip8(ss2 / ww);
                imOut->image[y][x*4 + 3] = clip8(ss3 / ww);
            }
        }
    }

    ImagingSectionLeave(&cookie);

    return imOut;
}

Imaging
ImagingTransform(Imaging imOut, Imaging imIn, int method,
                 int x0, int y0, int x1, int y1,
                 double a[8], int filterid, int fill)
{
    ImagingTransformMap transform;

    switch(method) {
    case IMAGING_TRANSFORM_AFFINE:
        return ImagingTransformAffineSlow(
            imOut, imIn, x0, y0, x1, y1, a, filterid, fill);
        break;
    case IMAGING_TRANSFORM_PERSPECTIVE:
        transform = perspective_transform;
        break;
    case IMAGING_TRANSFORM_QUAD:
        transform = quad_transform;
        break;
    default:
        return (Imaging) ImagingError_ValueError("bad transform method");
    }

    return ImagingGenericTransform(
        imOut, imIn,
        x0, y0, x1, y1,
        transform, a, filterid, fill);
}
