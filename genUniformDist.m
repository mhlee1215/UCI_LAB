function [ st ] = genUniformDist( ajkMat, e )
%GENUNIFORMDIST Summary of this function goes here
%   ai : 1 x T

    ajk_acc_t = zeros(size(ajkMat));
    for i=1:size(ajkMat, 1)
        ajk_acc_t(i,:) = sum(ajkMat(1:i, :), 1);
    end
    ajk_acc_s = zeros(size(ajkMat));
    ajk_acc_s(2:end, :) = ajk_acc_t(1:end-1, :);
    ajk_acc_tr = zeros(size(ajkMat));
    for i=1:size(ajkMat, 1)-1
        ajk_acc_tr(i, :) = sum(ajkMat(i+1:end, :), 1);
    end
    

    T = size(ajkMat, 1);
    K = size(ajkMat, 2);
    
    [s t] = meshgrid(1:T, 1:T);
    
    st = [];
    for i=1:K
        ajk = ajkMat(:,i);
        aj_acc_t = ajk_acc_t(:,i);
        aj_acc_tr = ajk_acc_tr(:,i);
        aj_acc_s = ajk_acc_s(:,i);
        
        aj_acc_tr2d = repmat(aj_acc_tr, 1, T); 
        aj_acc_t2d = repmat(aj_acc_t, 1, T);
        aj_acc_s2d = repmat(aj_acc_s', T, 1);
        
        %\sum_j\notin[s,t] \alpha_i log(\epsilon)
        p = (aj_acc_s2d.*log(e) + aj_acc_tr2d.*log(e)) .* (s <= t);
        %\sum_j\in[s,t] \alpha_i log(\gamma_{s,t})
%         real(log((1 - e.*(s - 1 - t + T))./(t-s+1)) .* (s <= t))
        p = p + real(log((1 - e.*(s - 1 - t + T))./(t-s+1)) .* (aj_acc_t2d - aj_acc_s2d) .* (s <= t));
        p(find(isnan(p))) = 0;
        p = p + (s > t).*-9999;
        [~, idx] = max(p(:));
        [mt, ms] = ind2sub([T, T], idx);
        st(end+1, :) = [ms mt];
    end
    
%     [s t] = meshgrid(1:T, 1:T);
%     dist = zeros(T, T);
%     dist(sub2ind([10 10], s, t)) = (s-1).*e + (T-t).*e + (1 - e.*(t - s - T +1)) ./ (t - s + 1);
%     dist(sub2ind([10 10], s, t)) = (1 - e.*(s - 1 - t + T));
%     dist(sub2ind([10 10], s, t)) = (s-1).*e + (T-t).*e;
%     
%     (s-1).*e + (T-t).*e + (1 - e.*(s - 1 - t + T))./(t-s+1) .* (t-s+1)
% 
%     
%     ps = [0;ajkMat(:,1)];
%     pt = [ajkMat(:,1);0];
%     p_max = 0;
%     for ss=1:10
%         for tt=1:10
%             if ss > tt 
%                 continue;
%             end
%             curP = sum(ps(1:ss))*(s).*e + sum(pt(tt+1:end))*(s).*e + sum(p
%         end
%     end

end

