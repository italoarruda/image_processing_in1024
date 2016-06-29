%Embassamento Manual de imagem sintetica
%imagem de entrada
I = imread('cameraman.tif');
subplot(1,2,1);imshow(I);title('Imagem Original');
r = 5;
c = 5;
h = zeros(r,c);
R = (r+1)/2;
C = (c+1)/2;
sigma = 1.25;

[X,Y] = meshgrid(-(r-1)/2:(r-1)/2, -(c-1)/2:(c-1)/2);
h = exp(-(X.^2 + Y.^2)/(2*sigma.^2));
h = h/sum(sum(h));

GaussBlur = imfilter(I,h);
subplot(1,2,2);imshow(GaussBlur);title('Imagem com Gaussian Blurred');

%Embassamento em Varios niveis
I = imread('cameraman.tif');

Iblur1 = imgaussfilt(I,2);
Iblur2 = imgaussfilt(I,4);
Iblur3 = imgaussfilt(I,8);
Iblur4 = imgaussfilt(I,16);

figure
imshow(I)
title('Imagem Original')

figure
imshow(Iblur1)
title('Imagem Suavizada, \sigma = 2')

figure
imshow(Iblur2)
title('Imagem Suavizada, \sigma = 4')

figure
imshow(Iblur3)
title('Imagem Suavizada, \sigma = 8')
imwrite(Iblur3, 'cameraman-base.tif','tif');
figure
imshow(Iblur4)
title('Imagem Suavizada, \sigma = 16')