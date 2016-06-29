%% Desembassar uma imagem usando Blind Deconvolution

rng default;

%Imagem de entrada
I = checkerboard(8);

%aplicação do filtro gaussiano
PSF = fspecial('gaussian',7,10);

V = .0001;
BlurredNoisy = imnoise(imfilter(I,PSF),'gaussian',0,V);

%%
% Cria um vetor de peso para especificar os pixels que estão inclusos no processamento
%%

WT = zeros(size(I));

WT(5:end-4,5:end-4) = 1;

INITPSF = ones(size(PSF));

%%
% Performa o blind deconvolution.
%%

[J P] = deconvblind(BlurredNoisy,INITPSF,20,10*sqrt(V),WT);

%%
% Apresenta os Resultados
%%

subplot(221);imshow(BlurredNoisy);
title('Imagem Borrada e com Ruido');

subplot(222);imshow(PSF,[]);
title('Verdadeiro PSF');

subplot(223);imshow(J);
title('Imagem Desembassada');

subplot(224);imshow(P,[]);
title('PSF Recuperado');