%  	Compute attenuation weighting
%   Input:
%       A:   image, size is [m, n]
%       alpha: Attenuation coefficient, 根据兰博比尔定律
%   Output:
%       W:   Weighting expresing depth-dependent attenuation
function [ Dw ] = attenuationWeighting2( img )
img = double(img);
img = medfilt2(img, [10, 10]);
gridentx = zeros(size(img));
gridenty1 = zeros(size(img));
gridenty2 = zeros(size(img));
gridentx(1:end-1, :) = img(1:end-1, :) - img(2:end, :);
gridentx(gridentx < 0) = 0;
% 对角线
gridenty1(2:end, 2:end) = img(1:end-1, 1:end-1) - img(2:end, 2:end);
gridenty1(gridenty1 < 0) = 0;
gridenty2(1:end-1, 1:end-1) = img(2:end, 2:end) - img(1:end-1, 1:end-1);
gridenty2(gridenty2 < 0) = 0;
grident = gridentx + 0.2*gridenty1 + 0.2*gridenty2;
%grident = 255 * grident ./ img;
Dw = cumsum(grident, 1);
Dw = ( Dw - min(Dw(:)) ) ./ ( max(Dw(:)) - min(Dw(:)) );    % 归一化范围到[0, 1]

end