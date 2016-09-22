function [R,t,X,S,a,pk,T,TAssigned, TXQ, vis, XC, XN, Xin] = joint_align( V,N,C, maxNumIter, clusterDensity )
%JOINT_ALIGN Summary of this function goes here
%   Detailed explanation goes here

updateVis = 0;
if nargin < 4
    maxNumIter = 100;
end

if nargin < 5
    clusterDensity = 2;
end

M = length(V);
% idx = transpose(1:M);
% strIdx = arrayfun(@(j) sprintf('view %d',j),idx,'uniformoutput',false);
% clrmap = {[1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0] ...
%     ; [1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
% clrmap = clrmap(1:M);

% markerSize = {};
% for ii=1:M
%     markerSize{end+1} = 1+ii;
% end
% markerSize = markerSize';
% markerSize = markerSize(1:M);

% marker_set = {'s', 'x', '.', '^'};
% marker = {};
% for i=1:M
%     marker{end+1} = marker_set{mod(i, 4)+1};
% end
% marker = marker';

% K = 300;
% randIdx = randperm(size(V{1}', 1));
% Xin = V{1}(:, randIdx(1:K));

va = [];
for i = 1:M
    va = [va V{i}];
end
[Xin,~,~,~,density]= uniformSubSample( va', clusterDensity);
Xin = Xin(:, find(density > 0));
K = size(Xin, 2);
% 
XC = [];
XN = [];
if isempty(C) && isempty(N) %ONly V
[R,t,X,S,a,pk,T,TAssigned, TXQ, vis, ~, ~] = jrmpc_soft_with_normal(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 1, 'updateVis', updateVis);
elseif isempty(C) %Only N and V
[R,t,X,S,a,pk,T,TAssigned, TXQ, vis, XC, ~] = jrmpc_soft_with_normal(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
    'normal', N, 'normalLambda', .1, 'updateVis', updateVis);    
    
elseif isempty(N) %C and V
    [R,t,X,S,a,pk,T,TAssigned, TXQ, vis, ~, XN] = jrmpc_soft_with_normal(V,Xin, ...
        'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
        'color', C, 'colorLambda', 0.1, 'updateVis', updateVis);
else
    [R,t,X,S,a,pk,T,TAssigned, TXQ, vis, XC, XN] = jrmpc_soft_with_normal(V,Xin, ...
    'maxNumIter',maxNumIter,'gamma',0.1, 'updatepriors', 0, ...
    'normal', N, 'normalLambda', .1, 'color', C, 'colorLambda', .1, 'updateVis', updateVis);    
end



end

