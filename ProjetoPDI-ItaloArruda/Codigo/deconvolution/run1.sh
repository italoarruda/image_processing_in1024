#!/bin/sh

images/cameraman.pgm psf/box9.pgm 0.001 -c .001
mv output.pgm cameraman_deconv.pgm

images/text.pgm psf/disc11.pgm 0.001 -t -e 5 -g
mv output.pgm text_deconv.pgm

images/face.ppm psf/disc7.pgm 0.001 -t -e 5 -g
mv output.ppm face_deconv.ppm

