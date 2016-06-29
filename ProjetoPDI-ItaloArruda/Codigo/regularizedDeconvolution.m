%Regularized deconvolution, deconvreg
I = im2double(imread('cameraman.tif'));
PSF = fspecial('gaussian',7,10);
V = .01;
BlurredNoisy = imfilter(I,PSF);

[J LAGRA] = deconvreg(BlurredNoisy,PSF);

subplot(121); imshow(BlurredNoisy);
title('Imagem Embassado');
subplot(122); imshow(J);
title('Imagem Desembassada');
imwrite(J, 'cameraman-regularDecon.tif','tif');