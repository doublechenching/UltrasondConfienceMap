clc;
clear all;
addpath './abus_matlab/code/dicm2nii';
volume = dicomread('H:\ConfidenceMapMatlab\941\20170503092639\US\1.2.276.0.7230010.3.1.4.1084411823.872.1493775799.49');
volume = squeeze(volume);
[height, width, depth] = size(volume);
threshold = 0.5;
alpha = 1.5; 
beta = 90;
gamma = 0.03;
tic
confidence1 = zeros(300, 1000, 1 , depth, 'uint8');   
for i = 1 : depth
    disp(i)
    img = volume(:, :, i);
    [map] = confMap(img, alpha, beta, gamma);
    dsc_map = dsc2(map);
    dsc_map = im2uint8(dsc_map); 
    confidence1(:, :, 1, i) = dsc_map;
end
confidence_volume = zeros(300, 1000, 1 , 1000, 'uint8');
for i = 1 : 300
    img = confidence1(i, :, 1, :);
    img = squeeze(img);
    img = imresize(img, [1000, 1000]);
    confidence_volume(i, :, 1, :) = img;
end
dsc_label = confidence_volume > threshold*255;
dsc_label = uint8(dsc_label);
tmp = squeeze(dsc_label);
re_label = imresize3(tmp, 0.25, 'nearest');
confidence_volume250 = squeeze(confidence_volume);
confidence_volume250 = imresize3(confidence_volume250, 0.25);
confidence_volume250 = reshape(confidence_volume250, [75, 250, 1, 250]);
re_label = reshape(re_label, [75, 250, 1, 250]);
toc
dicomwrite(confidence_volume250, 'dsc_confidence_941_10_250x75x250.DICOM');
dicomwrite(re_label, 'dsc_label_941_10_250x75x250.DICOM');
