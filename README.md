# Pillow-SIMD

Pillow-SIMD is "following" Pillow fork (which is PIL fork itself).

For more information about original Pillow, please
[read the documentation](http://pillow.readthedocs.io/),
[check the changelog](https://github.com/python-pillow/Pillow/blob/master/CHANGES.rst) and
[find out how to contribute](https://github.com/python-pillow/Pillow/blob/master/.github/CONTRIBUTING.md).


## Aims


## Status


## Performance

Hardware: Intel Core i5 4258U with AVX2.

Software: Ubuntu 14.04 64-bit.

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


## Contributing to Pillow-SIMD

**Important.** Pillow-SIMD and Pillow are two separate projects.
Please submit bugs and improvements not related to SIMD to 
[original Pillow](https://github.com/python-pillow/Pillow/issues/new).
All bugs and fixes in Pillow will appear in next Pillow-SIMD version automatically.
