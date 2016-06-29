%----------------------SSIM Experiments----------------------

ref = double(imread('20-cameraman-base.tif')); 
A = double(imread('21-cameraman-regularDecon.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(imread('22-cameraman-lucy.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(imread('23-cameraman-wiener.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(imread('24-cameraman-art.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);


ref = double(rgb2gray(imread('0-caracteres-base.png'))); 
A = double(rgb2gray(imread('1-caracteres-regularDecon.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(rgb2gray(imread('2-caracteres-lucy.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(rgb2gray(imread('3-caracteres-wiener.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(rgb2gray(imread('4-caracteres-art.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);


ref = double(rgb2gray(imread('10-girl-base.png'))); 
A = double(rgb2gray(imread('11-girl-regularDecon.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(rgb2gray(imread('12-girl-lucy.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(rgb2gray(imread('13-girl-wiener.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(rgb2gray(imread('14-girl-art.png')));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

ref = double(imread('30-lena-base.tif')); 
A = double(imread('31-lena-regularDecon.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(imread('32-lena-lucy.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(imread('33-lena-wiener.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

A = double(imread('34-lena-art.tif'));
[ssimval, ssimmap] = ssim(A,ref);
fprintf('The SSIM value is %0.4f.\n',ssimval);

%----------------------PSNR Experiments----------------------

ref = imread('20-cameraman-base.tif');

A1 = imread('21-cameraman-regularDecon.tif');
a1 = psnr(A1, ref);
B1 = imread('22-cameraman-lucy.tif');
a2 = psnr(B1, ref);
C1 = imread('23-cameraman-wiener.tif');
a3 = psnr(C1, ref);
D1 = imread('24-cameraman-art.tif');
a4 = psnr(D1, ref);

fprintf('The PSNR value for a1 is %0.4f.\n',a1);
fprintf('The PSNR value for a2 is %0.4f.\n',a2);
fprintf('The PSNR value for a3 is %0.4f.\n',a3);
fprintf('The PSNR value for a4 is %0.4f.\n',a4);

ref = imread('0-caracteres-base.png');

A2 = imread('1-caracteres-regularDecon.png');
b1 = psnr(A2, ref);
B2 = imread('2-caracteres-lucy.png');
b2 = psnr(B2, ref);
C2 = imread('3-caracteres-wiener.png');
b3 = psnr(C2, ref);
D2 = imread('4-caracteres-art.png');
b4 = psnr(D2, ref);

fprintf('The PSNR value for b1 is %0.4f.\n',b1);
fprintf('The PSNR value for b2 is %0.4f.\n',b2);
fprintf('The PSNR value for b3 is %0.4f.\n',b3);
fprintf('The PSNR value for b4 is %0.4f.\n',b4);


ref = imread('10-girl-base.png');

A3 = imread('11-girl-regularDecon.png');
c1 = psnr(A3, ref);
B3 = imread('12-girl-lucy.png');
c2 = psnr(B3, ref);
C3 = imread('13-girl-wiener.png');
c3 = psnr(C3, ref);
D3 = imread('14-girl-art.png');
c4 = psnr(D3, ref);

fprintf('The PSNR value for c1 is %0.4f.\n',c1);
fprintf('The PSNR value for c2 is %0.4f.\n',c2);
fprintf('The PSNR value for c3 is %0.4f.\n',c3);
fprintf('The PSNR value for c4 is %0.4f.\n',c4);

ref = imread('30-lena-base.tif');

A4 = imread('31-lena-regularDecon.tif');
d1 = psnr(A4, ref);
B4 = imread('32-lena-lucy.tif');
d2 = psnr(B4, ref);
C4 = imread('33-lena-wiener.tif');
d3 = psnr(C4, ref);
D4 = imread('34-lena-art.tif');
d4 = psnr(D4, ref);

fprintf('The PSNR value for d1 is %0.4f.\n',d1);
fprintf('The PSNR value for d2 is %0.4f.\n',d2);
fprintf('The PSNR value for d3 is %0.4f.\n',d3);
fprintf('The PSNR value for d4 is %0.4f.\n',d4);