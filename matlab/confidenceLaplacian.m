%  	Compute 6-Connected Laplacian for confidence estimation problem
%   六连接， 四连击加对角线方向， 应该是八邻域
%   Input:
%       P: 距离加权图像的坐标, 边界经过padding 0
%       A: 距离加权图像, 边界经过padding 0
%       beta:   Random walks parameter
%       gamma:  Horizontal penalty factor
%   Output:
%       map:    Confidence map for data
function L = confidenceLaplacian( P, A, beta, gamma )
[m, n] = size(P);
p = find(P);    % 大于0的坐标
i = P(p);       % 坐标向量, 非padding区域
j = P(p);       % Index vector
s = zeros(size(p)); % Entries vector, initially for diagonal和图像大小一致
% Vertical edges
for k = [-1 1]
   Q = P(p + k);    % 上下方坐标   
   q = find(Q);     % 去除边界
   ii = P( p(q) );  
   i = [i; ii];
   jj = Q(q);
   j = [j; jj];
   W = abs(A(p(ii)) - A(p(jj)));     % Intensity derived weight基于亮度的权重
   s = [s; W];
end

vl = numel(s);                      % 元素个数
% Diagonal edges    对角线和左右边界
for k = [(m-1) (m+1) (-m-1) (-m+1)]   
   Q = P(p+k);
   q = find(Q);    
   ii = P(p(q));
   i = [i; ii]; 
   jj = Q(q);
   j = [j; jj];
   W = abs(A(p(ii)) - A(p(jj))); % Intensity derived weight
   s = [s; W];
   
end


% Horizontal edges
for k = [m -m]
   Q = P(p+k);  
   q = find(Q);   
   ii = P(p(q));
   i = [i; ii];
   jj = Q(q);
   j = [j; jj];
   W = abs(A(p(ii))-A(p(jj))); % Intensity derived weight
   s = [s; W];
end
% s为亮度插值矩阵
% Normalize weights
s = (s - min(s(:))) ./ (max(s(:)) - min(s(:)) + eps);

% Horizontal penalty
s(vl+1:end) = s(vl+1:end) + gamma;

% Normalize differences
s = (s - min(s(:))) ./ (max(s(:)) - min(s(:)) + eps);

% Gaussian weighting function
EPSILON = 10e-6;
s = -((exp(-beta * s)) + EPSILON);
% Create Laplacian, diagonal missing
L = sparse(i, j, s); % i,j indices, s entry (non-zero) vector 创建稀疏矩阵L

% Reset diagonal weights to zero for summing 
% up the weighted edge degree in the next step, 因为对角元素，为结点本身，其值为边界权重之和
L = spdiags(zeros(size(s, 1), 1), 0, L);        % replaces the diagonals specified by 0 with the columns of zeros(size(s, 1). The output is sparse.
% Weighted edge degree
D = full(abs(sum(L, 1)))';                      % 结点值

% Finalize Laplacian by completing the diagonal
L = spdiags(D,0,L);                             % 将D复制回对角线
end

