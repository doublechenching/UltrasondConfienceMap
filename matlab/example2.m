clc;
clear all;
volume = dicomread('1.2.276.0.7230010.3.1.4.1719946590.156.1482464216.478');
volume = squeeze(volume);
img = volume(:, :, 150);
%img = imread('1.bmp');
threshold = 0.45;
alpha = 1.5;
beta = 120;
%gamma = 0.03;
gamma = 0.07;
[map] = confMap(img, alpha, beta, gamma);
label = map > threshold;
contour = bwperim(label);
figure;
subplot(311);
img(contour) = 255;
imshow(img);
subplot(312);
imshow(map);
subplot(313);
imshow(label);



