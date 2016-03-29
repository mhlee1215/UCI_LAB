function [ vis ] = genVisFromPeriod( st, dSize, K, e )
%GENVISFROMPERIOD Summary of this function goes here
%   Detailed explanation goes here
   visMat = zeros(dSize, K);
   T = dSize;
   
       for j=1:K
            %visMat(:,j) = 0.1;
            
            s = st(j,1);
            t = st(j,2);
            visMat(s:t,j) = (1 - e.*(s - 1 - t + T))./(t-s+1);
            visMat(find(visMat(:,j) == 0),j) = e;
       end
   
   
   vis = mat2cell(visMat, repmat(1, 1, dSize), K)';
end

