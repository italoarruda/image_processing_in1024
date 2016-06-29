%Lucy Richardson deconvolution, deconvlucy
I = im2double(imread('cameraman.tif'));
PSF = fspecial('gaussian',7,10);
V = .0001;
BlurredNoisy = imfilter(I,PSF);

J1 = deconvlucy(BlurredNoisy,PSF);

subplot(121);imshow(BlurredNoisy);
title('A = Embassado');
subplot(122);imshow(J1);
title('Deconvolucao Richard-Lucy');
imwrite(J1, 'cameraman-lucy.tif','tif');