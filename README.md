# Pillow-SIMD

Pillow-SIMD is "following" Pillow fork (which is PIL fork itself).

For more information about original Pillow, please
[read the documentation](http://pillow.readthedocs.io/),
[check the changelog](https://github.com/python-pillow/Pillow/blob/master/CHANGES.rst) and
[find out how to contribute](https://github.com/python-pillow/Pillow/blob/master/.github/CONTRIBUTING.md).


## What is SIMD


## Status


## Benchmarks

Numbers are megapixels of source image per second. For example if resize of
7712x4352 image is done in 0.5 seconds, the resoult will be 61.7 Mpx/s.

- ImageMagick 6.9.3-8 Q8 x86_64
- Pillow 3.2.0
- Pillow-SIMD 3.2.0.post1

Source    | Operation                   | IM   | Pillow | SIMD SSE4 | SIMD AVX2 
----------|-----------------------------|------|--------|-----------|-----------
7712x4352 | Bilinear resize to 16x16    | 27.0 | 217    | 456       | 545
          | Bilinear resize to 320x180  | 32.0 | 166    | 354       | 410
          | Bilinear resize to 2048x1155| 20.7 | 87.6   | 202       | 217
          | Bicubic resize to 16x16     | 10.9 | 115    | 240       | 278
          | Bicubic resize to 320x180   | 16.5 | 92.3   | 198       | 204
          | Bicubic resize to 2048x1155 | 12.2 | 65.7   | 126       | 130
          | Lanczos resize to 16x16     | 6.6  | 76.1   | 162       | 194
          | Lanczos resize to 320x180   | 11.0 | 63.2   | 133       | 147
          | Lanczos resize to 2048x1155 | 8.7  | 41.3   | 88.2      | 95.6
          | Blur 1px                    | 8.1  | 17.1   | 37.8
          | Blur 10px                   | 2.6  | 17.4   | 39.0
          | Blur 100px                  | 0.3  | 17.2   | 39.0
1920x1280 | Bilinear resize to 16x16    | 41.6 | 196    | 422       | 489
          | Bilinear resize to 320x180  | 27.6 | 111    | 196       | 197
          | Bilinear resize to 2048x1155| 9.1  | 20.7   | 71.3      | 72.6
          | Bicubic resize to 16x16     | 18.9 | 102    | 225       | 263
          | Bicubic resize to 320x180   | 14.5 | 66.3   | 154       | 162
          | Bicubic resize to 2048x1155 | 6.3  | 16.9   | 49.3      | 54.3
          | Lanczos resize to 16x16     | 13.7 | 68.6   | 118       | 167
          | Lanczos resize to 320x180   | 9.8  | 44.3   | 102       | 107
          | Lanczos resize to 2048x1155 | 4.7  | 14.6   | 36.8      | 40.6
          | Blur 1px                    | 8.7  | 16.2   | 35.7
          | Blur 10px                   | 2.8  | 16.7   | 35.4
          | Blur 100px                  | 0.4  | 16.4   | 36.2


### Methodology

All tests done on Ubuntu 14.04 64-bit runing on
Intel Core i5 4258U with AVX2 CPU.

Following script were using for tests:
https://gist.github.com/homm/f9b8d8a84a57a7e51f9c2a5828e40e63


## Why Pillow is so fast itself

There are no cheats. High-quality resize and blur methods are used for all
benchmakrs. And result is 


## Why Pillow-SIMD even faster


## Contributing to Pillow-SIMD

**Important.** Pillow-SIMD and Pillow are two separate projects.
Please submit bugs and improvements not related to SIMD to 
[original Pillow](https://github.com/python-pillow/Pillow/issues/new).
All bugs and fixes in Pillow will appear in next Pillow-SIMD version automatically.
