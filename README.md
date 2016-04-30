# Pillow-SIMD

Pillow-SIMD is "following" Pillow fork (which is PIL fork itself).

For more information about original Pillow, please
[read the documentation][original-docs],
[check the changelog][original-changelog] and
[find out how to contribute][original-contribute].


## Why SIMD

There are many ways to improve performance of image processing.
You can use better algorithms for the same task, you can make better
implementation for current algorithms, or you can use more processing unit
resources. It is perfect, when you can just use more efficient algirithm like
when gaussian blur based on convolutions [was replaced][gaussian-blur-changes]
with sequential box filters. But the number of such improvements is very limited.
It is also very tempting to use more processor unit resources 
(via parallelization), when they are available. But it is more handy just
to make things faster on the same resources. And that is where SIMD works better.

SIMD stands for "single instruction, multiple data". This is a way to perform
same operations against the huge amount of homogeneous data. 
Modern CPU have differnt SIMD instructions sets, like
MMX, SSE-SSE4, AVX, AVX2, AVX512, NEON.

Currently, Pillow-SIMD can be [compiled](#installation) with SSE4 (default)
and AVX2 support.


## Status

Pillow-SIMD can be used in production. Pillow-SIMD operates on
[Uploadcare](https://uploadcare.com/) servers for more than 1 year.
Uploadcare is SAAS for image storing and processing in the cloud
and the main sponsor of Pillow-SIMD project.

Cuurently, following methods are accelerated:

- Resize (convolustion-based resample): SSE4, AVX2
- Gaussian and box blur: SSE4


## Benchmarks

The numbers are processed megapixels of source image per second.
For example if resize of 7712x4352 image is done in 0.5 seconds,
the result will be 61.7 Mpx/s.

- ImageMagick 6.9.3-8 Q8 x86_64
- Pillow 3.2.0
- Pillow-SIMD 3.2.0.post1

Source    | Operation               | Filter  | IM   | Pillow | SIMD SSE4 | SIMD AVX2 
----------|-------------------------|---------|------|--------|-----------|-----------
7712x4352 | **Resize to 16x16**     | Bilinear| 27.0 | 217    | 456       | 545
          |                         | Bicubic | 10.9 | 115    | 240       | 278
          |                         | Lanczos | 6.6  | 76.1   | 162       | 194
          | **Resize to 320x180**   | Bilinear| 32.0 | 166    | 354       | 410
          |                         | Bicubic | 16.5 | 92.3   | 198       | 204
          |                         | Lanczos | 11.0 | 63.2   | 133       | 147
          | **Resize to 2048x1155** | Bilinear| 20.7 | 87.6   | 202       | 217
          |                         | Bicubic | 12.2 | 65.7   | 126       | 130
          |                         | Lanczos | 8.7  | 41.3   | 88.2      | 95.6
          | **Blur**                | 1px     | 8.1  | 17.1   | 37.8
          |                         | 10px    | 2.6  | 17.4   | 39.0
          |                         | 100px   | 0.3  | 17.2   | 39.0
1920x1280 | **Resize to 16x16**     | Bilinear| 41.6 | 196    | 422       | 489
          |                         | Bicubic | 18.9 | 102    | 225       | 263
          |                         | Lanczos | 13.7 | 68.6   | 118       | 167
          | **Resize to 320x180**   | Bilinear| 27.6 | 111    | 196       | 197
          |                         | Bicubic | 14.5 | 66.3   | 154       | 162
          |                         | Lanczos | 9.8  | 44.3   | 102       | 107
          | **Resize to 2048x1155** | Bilinear| 9.1  | 20.7   | 71.3      | 72.6
          |                         | Bicubic | 6.3  | 16.9   | 49.3      | 54.3
          |                         | Lanczos | 4.7  | 14.6   | 36.8      | 40.6
          | **Blur**                | 1px     | 8.7  | 16.2   | 35.7
          |                         | 10px    | 2.8  | 16.7   | 35.4
          |                         | 100px   | 0.4  | 16.4   | 36.2


### Some conclusion

Pillow is always faster than ImageMagick. And Pillow-SIMD is always faster
than Pillow in 2—2.5 time. In general Pillow-SIMD with AVX2 almost always
**10 times faster** than ImageMagick.

### Methodology

All tests done on Ubuntu 14.04 64-bit runing on
Intel Core i5 4258U with AVX2 CPU on single thread.

ImageMagick performance measured with command-line tool `convert` with
`-verbose` and `-bench` arguments. I'm using command line because
I needed to test latest version and this is to easist way to do that.
All operations are producing exactly the same result.
Resizing filters compliance:

- PIL.Image.BILINEAR == Triangle
- PIL.Image.BICUBIC == Catrom
- PIL.Image.LANCZOS == Lanczos

In ImageMagick the radius of gaussian blur called sigma and second parameter
is called radius. In fact there should not be additional parameters for
*gaussian blur*, because if the radius is too small, this is *not*
gaussian blur anymore. And if the radius is big, this does not give any
advantages, but makes operation slower. For test I set radius to sigma × 2.5.

Following script were used for testing:
https://gist.github.com/homm/f9b8d8a84a57a7e51f9c2a5828e40e63


## Why Pillow itself is so fast

There are no cheats. High-quality resize and blur methods are used for all
benchmakrs. And result is 


## Why Pillow-SIMD is even faster


## Why not to contribute SIMD to the original Pillow

Well, this is not so simple. First of all, Pillow support large number
of architectures, not only x86. But even for x86 platforms Pillow orfen
distributed via precompiled binaries. To integrate SIMD in precompiled binaries
we need to do runtime checks of CPU capabilites.
To compile code with runtime checks we need to pass `-mavx2` option
to compiller. But this automaticaly activates all `if (__AVX2__)` and below
conditions. And SIMD instructions under such conditions are exist
even in standard C library and they do not have any runtime checks.
Currently I don't know a way how to tell to compiller to allow SIMD
instrustions in the code but *do not allow* using such instructions without
runtime checks.


## Installation

In general you need to do `pip install pillow-simd` as always and if you
are using SSE4-capable CPU everything should run smoothly.
Do not forget to remove original Pillow package first.

If you want AVX2-enabled version, you need to pass additional flag to C
compiller. The easiest way to do that is define `CC` variable while compilation.

```bash
$ pip uninstall pillow
$ CC="cc -mavx2" pip install -U --force-reinstall pillow-simd
```


## Contributing to Pillow-SIMD

**Important.** Pillow-SIMD and Pillow are two separate projects.
Please submit bugs and improvements not related to SIMD to 
[original Pillow][original-issues]. All bugs and fixes in Pillow
will appear in next Pillow-SIMD version automatically.


  [original-docs]: http://pillow.readthedocs.io/
  [original-issues]: https://github.com/python-pillow/Pillow/issues/new
  [original-changelog]: https://github.com/python-pillow/Pillow/blob/master/CHANGES.rst
  [original-contribute]: https://github.com/python-pillow/Pillow/blob/master/.github/CONTRIBUTING.md
  [gaussian-blur-changes]: http://pillow.readthedocs.io/en/3.2.x/releasenotes/2.7.0.html#gaussian-blur-and-unsharp-mask
