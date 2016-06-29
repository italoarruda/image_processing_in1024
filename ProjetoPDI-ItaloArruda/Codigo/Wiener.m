%Wiener deconvolution, deconvwnr

I = im2double(imread('cameraman.tif'));
figure
imshow(I);
title('Imagem Original');
LEN = 21;
THETA = 11;
PSF = fspecial('motion', LEN, THETA);
blurred = imfilter(I, PSF, 'conv', 'circular');

figure
imshow(blurred);
title('Blurred Image');

figure
wnr1 = deconvwnr(blurred, PSF, 0);
imshow(wnr1);
title('Restored Image');
noise_mean = 0;
noise_var = 0.0001;
blurred_noisy = imnoise(blurred, 'gaussian', ...
                        noise_mean, noise_var);
imwrite(wnr1, 'cameraman-wiener.tif','tif');