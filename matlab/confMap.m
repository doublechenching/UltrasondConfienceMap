%  	Compute confidence map
%   Input:
%       data:   超声的RF信号数据，每一列为一个扫描线[m, n]
%       mode:   'RF' or 'B' mode data，模式选择，RF和B型超声数据
%       alpha, beta, gamma: See Medical Image Analysis reference根据图像类型调整的参数，见论文    
%   Output:
%       map:    Confidence map for data, 二维的矩阵的概率图
function [ map ] = confMap( data, alpha, beta, gamma, mode)
% 默认参数设置
if nargin < 4
    alpha = 2.0;
    beta = 90;
    gamma = 0.05;
end
% 默认图像类型
if(nargin < 5)
    mode = 'B';
end
% 图像归一化
data = double(data);
data = (data - min(data(:))) ./ ((max(data(:))-min(data(:)))+eps);
% 转换到希尔伯特空间v
if(strcmp(mode, 'RF'))
    data = abs(hilbert(data));
end
% Seeds and labels (boundary conditions)
seeds = [];
labels = [];
sc = 1: size(data, 2);                   %All elements，列数
sr_up = ones(1, length(sc));            % 注意这里全部为1
seed = sub2ind(size(data), sr_up, sc);  % 将下标转化为线性排列的坐标，这列转化多个。行坐标全部为1，代表所有第一行元素的坐标
seed = unique(seed);
seeds = [seeds seed];
% Label 1
label = zeros(1,length(seed));
label = label + 1;
labels = [labels label];                % 探头元素设置标签为1
sr_down = ones(1,length(sc));
sr_down = sr_down * size(data,1);       % 全部为[n_rows....n_rows]
seed = sub2ind(size(data),sr_down,sc);  % 所有最后一行的坐标
seed = unique(seed);
seeds = [seeds seed];
%Label 2
label = zeros(1,length(seed));
label = label + 2;
labels = [labels label];                % 最后一行标记为2
% Attenuation with Beer-Lambert
% W = attenuationWeighting(data, alpha);  % 衰减权重矩阵
W = attenuationWeighting2(data);
% Apply weighting directly to image
% Same as applying it individually during the formation of the Laplacian
data = data .* W;                       % 图像加权
% Find condidence values求解信心值
map = confidenceEstimation( data, seeds, labels, beta, gamma);
% Only keep probabilities for virtual source notes.
map = map(:, :, 1);   % size is [n_rows, n_cols, 2]
end

