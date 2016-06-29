Gradient domain deconvolution code

        from

Yosuke Bando and Tomoyuki Nishita
"Towards Digital Refocusing from a Single Photograph"
Proceedings of the 15th Pacific Conference on Computer Graphics
and Applications (Pacific Graphics 2007), pp. 363-372, 2007.

        and

United States Patent Application 20080175508.


                                  Apr. 24, 2011
                                  Yosuke Bando
                                  ybando@nis-lab.is.s.u-tokyo.ac.jp
                                  yosuke1.bando@toshiba.co.jp


0. Preface

The code is provided without any warranty, and it may be
used for non-commercial purposes for free.



1. To build

Makefile for Linux/Cygwin is included. Just type:

% make

An executable named "main(.exe)" will be created.



2. For test run

% ./run.sh

Three deconvolved images (.pgm/ppm) will be created.



3. Description of the files

+ images/         : sample blurred images
  + cameraman.pgm : well-known cameraman image
                    synthetically blurred with a 9x9 box kernel
  + face.ppm      : picture of a face
                    captured with Canon EOS-1D + EF 28-70mm lens
  + text.pgm      : picture of printed alphabets
                    captured with Canon EOS-1D + EF 28-70mm lens

+ psf/            : point-spread functions (PSF)
  + box9.pgm      : 9x9 box kernel used to blur the cameraman image
  + disc7.pgm     : disc kernel with a radius of 7 pixels
                    which approximates the defocus blur in face.ppm
  + disc11.pgm    : disc kernel with a radius of 11 pixels
                    which approximates the defocus blur in text.pgm

+ fft2d/          : 2D FFT code written by Takuya Ooura
                    (see fft2d/readme2d.txt for details)

+ *.{h,c}         : source code of the gradient domain deconvolution

+ Makefile        : Makefile for Linux/Cygwin

+ run.sh          : sample commands for running the executable

+ README.txt      : this file



4. Command line arguments

main in.pgm/ppm psf.pgm prior [options]
-e int  : edge taper size
-g      : undo gamma correction
-c float: convergence threshold
-n int  : maximum number of iterations
-t      : two-stage iteration
-r int  : post-process RL iterations

Mandatory arguments
(1) in.pgm/ppm: input blurred image (either in grayscale or color)
(2) psf.pgm   : point-spread function used to deconvolve the input image
(3) prior     : float value specifying the weight for the prior term,
                denoted by w in the paper, whose typical value is 0.001
                (larger values will better suppress ringing and noise
                but may also corrupt fine textures)

Optional arguments
-e int
    This option specifies the size of the edge taper in pixels.
    A region within the specified distance from the image border
    will be blurred by the PSF to reduce boundary effects of
    deconvolution. Default is 0. Values comparable to the size
    of a PSF are recommended.

-g
    In most cases, images captured with digital cameras are gamma-
    corrected, which needs to be inverted before deconvolution.
    The inversion is performend when this option is specified.

-c float
    This option specifies the convergence threshold for the iterative
    deconvolution. Default is 0.01. Smaller values may be specified for
    better results, although that will lead to longer execution time
    and could have adverse effects unless an accurate PSF is provided.

-n int
    This option specifies the maximum number of deconvolution iterations.
    Default is infinity. Often, two or three iterations are sufficient.

-t
    Better results may be obtained by first setting a larger prior weight
    to suppress ringing, and then decreasing it to put more weight on
    the data fidelity term. When this option is specified, deconvolution
    is first performed with the prior weight that is twice as large as
    the specified value, and then with the original value.

-r int
    Gradient domain deconvolution does not enforce positivity constraint
    on pixel values. Possible ramification of this, if any, may be
    ameliorated by some post-process Richardson-Lucy (RL) iterations.
    This option specifies the number of RL iterations. Default is 0.

