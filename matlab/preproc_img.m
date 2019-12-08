function [img] = preproc_img(path)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
img = imread(path);
img = im2double(img);
img = imadjust(img, [0.007, 0.81]);

end

