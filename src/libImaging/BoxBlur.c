#include "Python.h"
#include "Imaging.h"


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


typedef UINT8 pixel[4];

/* General implementation when radius > 0 and not too big */
void static inline
ImagingLineBoxBlur4(pixel *lineOut, pixel *lineIn, int lastx, int radius,
    int edgeA, int edgeB, UINT32 ww, UINT32 fw)
{
    int x;
    UINT32 acc[4];
    UINT32 bulk[4];

    #define MOVE_ACC(acc, subtract, add) \
        acc[0] += lineIn[add][0] - lineIn[subtract][0]; \
        acc[1] += lineIn[add][1] - lineIn[subtract][1]; \
        acc[2] += lineIn[add][2] - lineIn[subtract][2]; \
        acc[3] += lineIn[add][3] - lineIn[subtract][3];

    #define ADD_FAR(bulk, acc, left, right) \
        bulk[0] = (acc[0] * ww) + (lineIn[left][0] + lineIn[right][0]) * fw; \
        bulk[1] = (acc[1] * ww) + (lineIn[left][1] + lineIn[right][1]) * fw; \
        bulk[2] = (acc[2] * ww) + (lineIn[left][2] + lineIn[right][2]) * fw; \
        bulk[3] = (acc[3] * ww) + (lineIn[left][3] + lineIn[right][3]) * fw;

    #define SAVE(x, bulk) \
        lineOut[x][0] = (UINT8)((bulk[0] + (1 << 23)) >> 24); \
        lineOut[x][1] = (UINT8)((bulk[1] + (1 << 23)) >> 24); \
        lineOut[x][2] = (UINT8)((bulk[2] + (1 << 23)) >> 24); \
        lineOut[x][3] = (UINT8)((bulk[3] + (1 << 23)) >> 24);

    /* Compute acc for -1 pixel (outside of image):
       From "-radius-1" to "-1" get first pixel,
       then from "0" to "radius-1". */
    acc[0] = lineIn[0][0] * (radius + 1);
    acc[1] = lineIn[0][1] * (radius + 1);
    acc[2] = lineIn[0][2] * (radius + 1);
    acc[3] = lineIn[0][3] * (radius + 1);
    /* As radius can be bigger than xsize, iterate to edgeA -1. */
    for (x = 0; x < edgeA - 1; x++) {
        acc[0] += lineIn[x][0];
        acc[1] += lineIn[x][1];
        acc[2] += lineIn[x][2];
        acc[3] += lineIn[x][3];
    }
    /* Then multiply remainder to last x. */
    acc[0] += lineIn[lastx][0] * (radius - edgeA + 1);
    acc[1] += lineIn[lastx][1] * (radius - edgeA + 1);
    acc[2] += lineIn[lastx][2] * (radius - edgeA + 1);
    acc[3] += lineIn[lastx][3] * (radius - edgeA + 1);

    if (edgeA <= edgeB)
    {
        /* Subtract pixel from left ("0").
           Add pixels from radius. */
        for (x = 0; x < edgeA; x++) {
            MOVE_ACC(acc, 0, x + radius);
            ADD_FAR(bulk, acc, 0, x + radius + 1);
            SAVE(x, bulk);
        }
        /* Subtract previous pixel from "-radius".
           Add pixels from radius. */
        for (x = edgeA; x < edgeB; x++) {
            MOVE_ACC(acc, x - radius - 1, x + radius);
            ADD_FAR(bulk, acc, x - radius - 1, x + radius + 1);
            SAVE(x, bulk);
        }
        /* Subtract previous pixel from "-radius".
           Add last pixel. */
        for (x = edgeB; x <= lastx; x++) {
            MOVE_ACC(acc, x - radius - 1, lastx);
            ADD_FAR(bulk, acc, x - radius - 1, lastx);
            SAVE(x, bulk);
        }
    }
    else
    {
        for (x = 0; x < edgeB; x++) {
            MOVE_ACC(acc, 0, x + radius);
            ADD_FAR(bulk, acc, 0, x + radius + 1);
            SAVE(x, bulk);
        }
        for (x = edgeB; x < edgeA; x++) {
            MOVE_ACC(acc, 0, lastx);
            ADD_FAR(bulk, acc, 0, lastx);
            SAVE(x, bulk);
        }
        for (x = edgeA; x <= lastx; x++) {
            MOVE_ACC(acc, x - radius - 1, lastx);
            ADD_FAR(bulk, acc, x - radius - 1, lastx);
            SAVE(x, bulk);
        }
    }

    #undef MOVE_ACC
    #undef ADD_FAR
    #undef SAVE
}


/* Optimized implementation when radius = 0 */
void static inline
ImagingLineBoxBlur4Zero(pixel *lineOut, pixel *lineIn, int lastx,
    int edgeA, int edgeB, UINT32 ww, UINT32 fw)
{
    int x;
    UINT32 acc[4];
    UINT32 bulk[4];

    #define MOVE_ACC(acc, add) \
        acc[0] = lineIn[add][0]; \
        acc[1] = lineIn[add][1]; \
        acc[2] = lineIn[add][2]; \
        acc[3] = lineIn[add][3];

    #define ADD_FAR(bulk, acc, left, right) \
        bulk[0] = (acc[0] * ww) + (lineIn[left][0] + lineIn[right][0]) * fw; \
        bulk[1] = (acc[1] * ww) + (lineIn[left][1] + lineIn[right][1]) * fw; \
        bulk[2] = (acc[2] * ww) + (lineIn[left][2] + lineIn[right][2]) * fw; \
        bulk[3] = (acc[3] * ww) + (lineIn[left][3] + lineIn[right][3]) * fw;

    #define SAVE(x, bulk) \
        lineOut[x][0] = (UINT8)((bulk[0] + (1 << 23)) >> 24); \
        lineOut[x][1] = (UINT8)((bulk[1] + (1 << 23)) >> 24); \
        lineOut[x][2] = (UINT8)((bulk[2] + (1 << 23)) >> 24); \
        lineOut[x][3] = (UINT8)((bulk[3] + (1 << 23)) >> 24);

    if (edgeA <= edgeB)
    {
        for (x = 0; x < edgeA; x++) {
            MOVE_ACC(acc, x);
            ADD_FAR(bulk, acc, 0, x + 1);
            SAVE(x, bulk);
        }
        for (x = edgeA; x < edgeB; x++) {
            MOVE_ACC(acc, x);
            ADD_FAR(bulk, acc, x - 1, x + 1);
            SAVE(x, bulk);
        }
        for (x = edgeB; x <= lastx; x++) {
            MOVE_ACC(acc, lastx);
            ADD_FAR(bulk, acc, x - 1, lastx);
            SAVE(x, bulk);
        }
    }
    else
    {
        /* This is possible when radius = 0 and width = 1. */
        for (x = edgeB; x < edgeA; x++) {
            MOVE_ACC(acc, lastx);
            ADD_FAR(bulk, acc, 0, lastx);
            SAVE(x, bulk);
        }
    }

    #undef MOVE_ACC
    #undef ADD_FAR
    #undef SAVE
}


/* General implementation when radius > 0 and not too big */
void static inline
ImagingLineBoxBlur1(UINT8 *lineOut, UINT8 *lineIn, int lastx, int radius,
                    int edgeA, int edgeB, UINT32 ww, UINT32 fw)
{
    int x;
    UINT32 acc;
    UINT32 bulk;

    #define MOVE_ACC(acc, subtract, add) \
        acc += lineIn[add] - lineIn[subtract];

    #define ADD_FAR(bulk, acc, left, right) \
        bulk = (acc * ww) + (lineIn[left] + lineIn[right]) * fw;

    #define SAVE(x, bulk) \
        lineOut[x] = (UINT8)((bulk + (1 << 23)) >> 24)

    acc = lineIn[0] * (radius + 1);
    for (x = 0; x < edgeA - 1; x++) {
        acc += lineIn[x];
    }
    acc += lineIn[lastx] * (radius - edgeA + 1);

    if (edgeA <= edgeB)
    {
        for (x = 0; x < edgeA; x++) {
            MOVE_ACC(acc, 0, x + radius);
            ADD_FAR(bulk, acc, 0, x + radius + 1);
            SAVE(x, bulk);
        }
        for (x = edgeA; x < edgeB; x++) {
            MOVE_ACC(acc, x - radius - 1, x + radius);
            ADD_FAR(bulk, acc, x - radius - 1, x + radius + 1);
            SAVE(x, bulk);
        }
        for (x = edgeB; x <= lastx; x++) {
            MOVE_ACC(acc, x - radius - 1, lastx);
            ADD_FAR(bulk, acc, x - radius - 1, lastx);
            SAVE(x, bulk);
        }
    }
    else
    {
        for (x = 0; x < edgeB; x++) {
            MOVE_ACC(acc, 0, x + radius);
            ADD_FAR(bulk, acc, 0, x + radius + 1);
            SAVE(x, bulk);
        }
        for (x = edgeB; x < edgeA; x++) {
            MOVE_ACC(acc, 0, lastx);
            ADD_FAR(bulk, acc, 0, lastx);
            SAVE(x, bulk);
        }
        for (x = edgeA; x <= lastx; x++) {
            MOVE_ACC(acc, x - radius - 1, lastx);
            ADD_FAR(bulk, acc, x - radius - 1, lastx);
            SAVE(x, bulk);
        }
    }

    #undef MOVE_ACC
    #undef ADD_FAR
    #undef SAVE
}


/* Optimized implementation when radius = 0 */
void static inline
ImagingLineBoxBlur1Zero(UINT8 *lineOut, UINT8 *lineIn, int lastx,
                        int edgeA, int edgeB, UINT32 ww, UINT32 fw)
{
    int x;
    UINT32 acc;
    UINT32 bulk;

    #define MOVE_ACC(acc, add) \
        acc = lineIn[add];

    #define ADD_FAR(bulk, acc, left, right) \
        bulk = (acc * ww) + (lineIn[left] + lineIn[right]) * fw;

    #define SAVE(x, bulk) \
        lineOut[x] = (UINT8)((bulk + (1 << 23)) >> 24)

    if (edgeA <= edgeB) {
        for (x = 0; x < edgeA; x++) {
            MOVE_ACC(acc, x);
            ADD_FAR(bulk, acc, 0, x + 1);
            SAVE(x, bulk);
        }
        for (x = edgeA; x < edgeB; x++) {
            MOVE_ACC(acc, x);
            ADD_FAR(bulk, acc, x - 1, x + 1);
            SAVE(x, bulk);
        }
        for (x = edgeB; x <= lastx; x++) {
            MOVE_ACC(acc, lastx);
            ADD_FAR(bulk, acc, x - 1, lastx);
            SAVE(x, bulk);
        }
    } else {
        /* This is possible when radius = 0 and width = 1. */
        for (x = edgeB; x < edgeA; x++) {
            MOVE_ACC(acc, lastx);
            ADD_FAR(bulk, acc, 0, lastx);
            SAVE(x, bulk);
        }
    }

    #undef MOVE_ACC
    #undef ADD_FAR
    #undef SAVE
}


Imaging
ImagingHorizontalBoxBlur(Imaging imOut, Imaging imIn, float floatRadius)
{
    ImagingSectionCookie cookie;

    int y;

    int radius = (int) floatRadius;
    UINT32 ww = (UINT32) (1 << 24) / (floatRadius * 2 + 1);
    UINT32 fw = ((1 << 24) - (radius * 2 + 1) * ww) / 2;

    int edgeA = MIN(radius + 1, imIn->xsize);
    int edgeB = MAX(imIn->xsize - radius - 1, 0);

    // printf(">>> %d, from %d to %d. %d %d\n", radius, edgeA, edgeB, ww, fw);

    ImagingSectionEnter(&cookie);

    if (imIn->image8)
    {
        for (y = 0; y < imIn->ysize; y++) {
            if (radius == 0) {
                ImagingLineBoxBlur1Zero(
                    imOut->image8[y], imIn->image8[y],
                    imIn->xsize - 1, edgeA, edgeB, ww, fw);
            } else {
                ImagingLineBoxBlur1(
                    imOut->image8[y], imIn->image8[y],
                    imIn->xsize - 1, radius, edgeA, edgeB, ww, fw);
            }
        }
    }
    else
    {
        for (y = 0; y < imIn->ysize; y++) {
            if (radius == 0) {
                ImagingLineBoxBlur4Zero(
                    (pixel *) imOut->image32[y], (pixel *) imIn->image32[y],
                    imIn->xsize - 1, edgeA, edgeB, ww, fw);
            } else {
                ImagingLineBoxBlur4(
                    (pixel *) imOut->image32[y], (pixel *) imIn->image32[y],
                    imIn->xsize - 1, radius, edgeA, edgeB, ww, fw);
            }
        }
    }

    ImagingSectionLeave(&cookie);

    return imOut;
}


/* General implementation when radius > 0 and not too big */
void
ImagingInnerVertBoxBlur(Imaging imOut, Imaging imIn, int lasty, int radius,
                        int edgeA, int edgeB, UINT32 ww, UINT32 fw) 
{
    #define LINE 1024

    int x, xx, y;
    int line = LINE;
    UINT32 acc[LINE];
    
    UINT8 *lineOut, *lineAdd, *lineLeft, *lineRight;
    

    #define INNER_LOOP(line)  \
        for (x = 0; x < line; x++) { \
            UINT32 bulk; \
            acc[x] += lineAdd[xx+x] - lineLeft[xx+x]; \
            bulk = (acc[x] * ww) + (lineLeft[xx+x] + lineRight[xx+x]) * fw; \
            lineOut[xx+x] = (UINT8)((bulk + (1 << 23)) >> 24); \
        }

    for (xx = 0; xx < imIn->linesize; xx += LINE) {
        if (xx + LINE > imIn->linesize) {
            line = imIn->linesize - xx;
        }
        lineLeft = (UINT8 *)imIn->image[0];
        lineRight = (UINT8 *)imIn->image[lasty];
        for (x = 0; x < line; x++) {
            acc[x] = (lineLeft[xx+x] * (radius + 1) +
                      lineRight[xx+x] * (radius - edgeA + 1));
        }
        for (y = 0; y < edgeA - 1; y++) {
            lineAdd = (UINT8 *)imIn->image[y];
            for (x = 0; x < line; x++) {
                acc[x] += lineAdd[xx+x];
            }
        }

        if (edgeA <= edgeB) {
            for (y = 0; y < edgeA; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineAdd = (UINT8 *)imIn->image[y + radius];
                lineLeft = (UINT8 *)imIn->image[0];
                lineRight = (UINT8 *)imIn->image[y + radius + 1];
                INNER_LOOP(line);
            }
            for (y = edgeA; y < edgeB; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineAdd = (UINT8 *)imIn->image[y + radius];
                lineLeft = (UINT8 *)imIn->image[y - radius - 1];
                lineRight = (UINT8 *)imIn->image[y + radius + 1];
                INNER_LOOP(line);
            }
            for (y = edgeB; y <= lasty; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineLeft = (UINT8 *)imIn->image[y - radius - 1];
                lineAdd = lineRight = (UINT8 *)imIn->image[lasty];
                INNER_LOOP(line);
            }
        } else {
            for (y = 0; y < edgeB; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineAdd = (UINT8 *)imIn->image[y + radius];
                lineLeft = (UINT8 *)imIn->image[0];
                lineRight = (UINT8 *)imIn->image[y + radius + 1];
                INNER_LOOP(line);
            }
            for (y = edgeB; y < edgeA; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineLeft = (UINT8 *)imIn->image[0];
                lineAdd = lineRight = (UINT8 *)imIn->image[lasty];
                INNER_LOOP(line);
            }
            for (y = edgeA; y <= lasty; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineLeft = (UINT8 *)imIn->image[y - radius - 1];
                lineAdd = lineRight = (UINT8 *)imIn->image[lasty];
                INNER_LOOP(line);
            }
        }
    }

    #undef INNER_LOOP
    #undef LINE
}


/* Optimized implementation when radius = 0 */
void
ImagingInnerVertBoxBlurZero(Imaging imOut, Imaging imIn, int lasty,
                            int edgeA, int edgeB, UINT32 ww, UINT32 fw)
{
    #define LINE 1024
    
    int x, xx, y;
    int line = LINE;
    
    UINT8 *lineOut, *lineAdd, *lineLeft, *lineRight;
    

    #define INNER_LOOP(line)  \
        for (x = 0; x < line; x++) { \
            UINT32 bulk; \
            bulk = (lineAdd[xx+x] * ww) + (lineLeft[xx+x] + lineRight[xx+x]) * fw; \
            lineOut[xx+x] = (UINT8)((bulk + (1 << 23)) >> 24); \
        }

    for (xx = 0; xx < imIn->linesize; xx += LINE) {
        if (xx + LINE > imIn->linesize) {
            line = imIn->linesize - xx;
        }

        if (edgeA <= edgeB) {
            for (y = 0; y < edgeA; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineAdd = (UINT8 *)imIn->image[y];
                lineLeft = (UINT8 *)imIn->image[0];
                lineRight = (UINT8 *)imIn->image[y + 1];
                INNER_LOOP(line);
            }
            for (y = edgeA; y < edgeB; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineAdd = (UINT8 *)imIn->image[y];
                lineLeft = (UINT8 *)imIn->image[y - 1];
                lineRight = (UINT8 *)imIn->image[y + 1];
                INNER_LOOP(line);
            }
            for (y = edgeB; y <= lasty; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineLeft = (UINT8 *)imIn->image[y - 1];
                lineAdd = lineRight = (UINT8 *)imIn->image[lasty];
                INNER_LOOP(line);
            }
        } else {
            for (y = edgeB; y < edgeA; y++) {
                lineOut = (UINT8 *)imOut->image[y];
                lineLeft = (UINT8 *)imIn->image[0];
                lineAdd = lineRight = (UINT8 *)imIn->image[lasty];
                INNER_LOOP(line);
            }
        }
    }

    #undef INNER_LOOP
    #undef LINE
}


Imaging
ImagingVerticalBoxBlur(Imaging imOut, Imaging imIn, float floatRadius)
{
    ImagingSectionCookie cookie;

    int radius = (int) floatRadius;
    UINT32 ww = (UINT32) (1 << 24) / (floatRadius * 2 + 1);
    UINT32 fw = ((1 << 24) - (radius * 2 + 1) * ww) / 2;

    int edgeA = MIN(radius + 1, imIn->ysize);
    int edgeB = MAX(imIn->ysize - radius - 1, 0);
    int lasty = imIn->ysize - 1;

    // printf(">>> %d, from %d to %d. %d %d\n", radius, edgeA, edgeB, ww, fw);

    ImagingSectionEnter(&cookie);

    if (radius == 0) {
        ImagingInnerVertBoxBlurZero(imOut, imIn, lasty,
                                    edgeA, edgeB, ww, fw);
    } else {
        ImagingInnerVertBoxBlur(imOut, imIn, lasty, radius,
                                edgeA, edgeB, ww, fw);
    }

    ImagingSectionLeave(&cookie);

    return imOut;
}


Imaging
ImagingBoxBlur(Imaging imOut, Imaging imIn, float radius, int n)
{
    int i;
    /* Real allocated buffer image */
    Imaging imBuffer;
    /* Pointers to current working images */
    Imaging imFrom, imTo;
    /* Temp pointer for swapping current working images */
    Imaging imTemp;

    if (n < 1) {
        return ImagingError_ValueError(
            "number of passes must be greater than zero"
        );
    }

    if (strcmp(imIn->mode, imOut->mode) ||
        imIn->type  != imOut->type  ||
        imIn->bands != imOut->bands ||
        imIn->xsize != imOut->xsize ||
        imIn->ysize != imOut->ysize)
        return ImagingError_Mismatch();

    if (imIn->type != IMAGING_TYPE_UINT8)
        return ImagingError_ModeError();

    if (!(strcmp(imIn->mode, "RGB") == 0 ||
          strcmp(imIn->mode, "RGBA") == 0 ||
          strcmp(imIn->mode, "RGBa") == 0 ||
          strcmp(imIn->mode, "RGBX") == 0 ||
          strcmp(imIn->mode, "CMYK") == 0 ||
          strcmp(imIn->mode, "L") == 0 ||
          strcmp(imIn->mode, "LA") == 0 ||
          strcmp(imIn->mode, "La") == 0))
        return ImagingError_ModeError();

    imBuffer = ImagingNewDirty(imIn->mode, imIn->xsize, imIn->ysize);
    if ( ! imBuffer)
        return NULL;

    /* We always do 2n passes. Odd to imBuffer, even to imOut */
    imTo = imBuffer;
    imFrom = imOut;

    /* First pass, use imIn instead of imFrom. */
    ImagingVerticalBoxBlur(imTo, imIn, radius);
    for (i = 1; i < n; i ++) {
        /* Swap current working images */
        imTemp = imTo; imTo = imFrom; imFrom = imTemp;
        ImagingVerticalBoxBlur(imTo, imFrom, radius);
    }

    /* Reuse imOut as a source and destination there. */
    for (i = 0; i < n; i ++) {
        /* Swap current working images */
        imTemp = imTo; imTo = imFrom; imFrom = imTemp;
        ImagingHorizontalBoxBlur(imTo, imFrom, radius);
    }

    ImagingDelete(imBuffer);

    return imOut;
}


Imaging ImagingGaussianBlur(Imaging imOut, Imaging imIn, float radius,
    int passes)
{
    float sigma2, L, l, a;

    sigma2 = radius * radius / passes;
    // from http://www.mia.uni-saarland.de/Publications/gwosdek-ssvm11.pdf
    // [7] Box length.
    L = sqrt(12.0 * sigma2 + 1.0);
    // [11] Integer part of box radius.
    l = floor((L - 1.0) / 2.0);
    // [14], [Fig. 2] Fractional part of box radius.
    a = (2 * l + 1) * (l * (l + 1) - 3 * sigma2);
    a /= 6 * (sigma2 - (l + 1) * (l + 1));

    return ImagingBoxBlur(imOut, imIn, l + a, passes);
}
