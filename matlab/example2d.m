clc;
clear all;
img = imread('1.bmp');
img = im2double(img);
% img = imadjust(img, [0.007, 0.81]);
alpha = 1.5;
beta = 120;
gamma = 0.08;
threshold = 0.45;
tic
rimg = imresize(img, 0.25, 'nearest');
[map] = confMap(rimg, alpha, beta, gamma);
map = imresize(map, size(img));
toc
figure();
subplot(221)
imshow(img);
subplot(223)
imshow(map);
subplot(222)
mask = double(map > threshold);
imshow(map > threshold);
subplot(224)
imshow(mask.*img);
res = [map map>threshold mask.*img ];
% imwrite(res, 'res2.bmp')