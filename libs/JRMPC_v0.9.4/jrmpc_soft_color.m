function [R,t,X,Q,a,pk,T, TAssign, TXQ] = jrmpc_soft_color(V,X,varargin)
%         JRMPC    Joint Registration of Multiple Point Clouds.
%            [R,t] = JRMPC(V,X) estimates the Euclidean transformation parameters R,t in order to rigidly register the views in V.
%            V is an M x 1 cell array of views, with each view V{j} j = 1:M represented as a 3 x Nj matrix of (cartesian) coordinates.
%            More specifically a column V{j}(:,i) of an element V{j} in V contains the coordinates [x y z]^T of i-th 3D point
%            of the j-th view. X is a 3 x K matrix of initial centers.
%
%            The algorithm returns two, M x 1 cell arrays, R and t. The elements of R are 3 x 3 rotation matrices, and the elements
%            of t are the 3 x 1 translation vectors. The pair R{j},t{j} corresponds to the view stored in V{j}.
%
%            [R,t,X] = JRMPC(V,X) returns also a 3 x K matrix X with the final means of the GMM components.
%
%            [R,t,X,S] = JRMPC(V,X) returns also a K x 1 vector S containing the final (isotropic) variances sigma_k^2 of the
%            GMM components in X. Note that S is computed only if it is requested as the algorithm works with precisions.
%
%            [R,t,X,S,a] = JRMPC(V,X) returns also an M x 1 cell a with the posterior probabilities of assignment. More specifically,
%            each element a{j} in a, is an Nj x K matrix where each element a{j}(i,k) in a{j} denotes p(z_{ji} = k | v_{ji}) computed in
%            the last iteration. Outlier probabilites are not returned as they are never computed explicitly.
%            They can be found as 1-sum(a{j},2).
%
%            [R,t,X,S,a,p] = JRMPC(V,X) returns also a K x 1 vector p containing the (final) priors of the GMM components.
%            if updatePriors = 0, p is equal to the initial mixing coefficients - priors.
%
%            [R,t,X,S,a,p,T] = JRMPC(V,X) returns a M x 2 x maxNumIter cell array T with the history of transformation parameters.
%            In more detail, T stores the R and t cells computed at each iteration, e.g at iter-th iteration R_iter = T(:,1,iter) and
%            t_iter = T(:,2,iter). T is allocated only if it is requested.
%
%            JRMPC(V,X,PARAM1,VAL1,PARAM2,VAL2,...) calls jrmpc specifying optional parameters (see below)
%            and corresponding values that control initialization. Parameter names are case insensitive.
%            Any number of parameters can be provided, their order does not matter, although specifying a parameter multiple times,
%            even with the same value, results in an error.
%
%            Parameters include:
%
%            'R'                      M x 1 cell specifying initial rotation matrices for each view.
%                                     The elements of R must be 3 x 3 matrices, orthonormality is not verified internally.
%                                     Default value: All R{j} are initialized with the identity matrix if no stated otherwise.
%
%            't'                      M x 1 cell with 3 x 1 vectors used as initialization for the translation of each view.
%                                     Elements of R and t must be non empty when initialization is provided.
%                                     Default value: Each t{j} is initialized with the arithmetic mean of V{j}, e.g -mean(V{j},2).
%
%            'S'                      Initial variances (sigma_k^2) for the GMM components, can be a K x 1 vector or a scalar.
%                                     If scalar is provided then all K components are initialized with the same variance.
%                                     Default value: All variances are initialized with the same value,
%                                     which is computed as the squared length of the diagonal of the bounding box that
%                                     contains all points of V, after aplying initial rototranslation.
%
%            'maxNumIter'             Specifies the number of iterations, runned iterations are from 0 to floor(maxNumIter-1).
%                                     Default value: 100.
%
%            'epsilon'                Artificial covariance flatten. A positive number added to S, after its update, at every iteration.
%                                     Default value: 1e-6.
%
%            'initialPriors'          Specifies the prior probabilities p of the GMM components, and implicitly defines the prior p_{K+1}
%                                     for the outlier class. It can be a K x 1 vector or a scalar. If p is scalar then that same value is
%                                     used for all components. The sum of all elements in p (or K*p if p is scalar), must be less than 1
%                                     as they represent a probability mass. p_{K+1} is computed internally as 1 - sum(p) if p is a vector,
%                                     or as p_{K+1} = 1-K*p otherwise. gamma is uniquely defined from p_{K+1} as 1 = (gamma+1)*sum(p).
%                                     Default value: The distribution of p_k is initialized as a uniform as p_k = 1/(K+1), k=1:K+1.
%
%            'gamma'                  Positive scalar specifying the outlier proportion in V. Used to compute the prior probability
%                                     p_{K+1} of the outlier component as gamma*sum_k(p_k). If gamma is provided then pk's are
%                                     initialized uniformly as sum_k(p_k) = 1/(gamma+1) => p_k = 1/(K*(gamma+1)). Paramater gamma is a
%                                     shortcut to set initialPriors uniformly, and therefore, either  'gamma' or 'initialPriors'
%                                     should be given at a time.
%                                     Default value: 1/K.
%
%            'updatePriors'           It is a flag that controls the update of p across iterations. The algorithm expects a scalar.
%                                     If it is (numeric) 0 then p is kept fixed otherwise priors are updated at every iteration.
%                                     Default value: 0.
%
%             References:
%               [1] Georgios D. Evangelidis, D. Kounades-Bastian, R. Horaud, and E.Z Psarakis,
%                   A Generative Model for the Joint Registration of Multiple Point Sets, ECCV, 2014.
%
%            $Revision: 0.9.4   $  $DATE: 24/05/2015 13:45 PM $

sqe = @(Y,X) pdist2(X', Y')'.^2;% sum(bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1])).^2,3);
sqe2 = @(Y,X) bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1]));

% ======================================================================== %
%                           C H E C K S                                    %
% ======================================================================== %

V = V(:);
M = numel(V);

[dim,K] = size(X);

if dim ~= 3
    error('X must be a 3 x K matrix.');
end

if mod(numel(varargin),2)
    error('odd number of optional parameters, Opt parames must be given as string-value pairs.')
end

for j=1:M
    if size(V{j},1) ~= 3
        error('V must be an M x 1 cell of 3 x .. matrices, V{%d} has %d in dimension 1.',j,size(V{j},1));
    end
end

% ======================================================================== %
%                       V A R A R G I N  P A R S E R                       %
% ======================================================================== %

isSetR = 0;
isSetT = 0;
isSetQ = 0;
isSetMaxNumIter = 0;
isSetInitialPriorsOrGamma = 0;
isSetEpsilon = 0;
isSetEpsilonST = 0;
isSetUpdatePriors = 0;
isSetVis = 0;
isSetUpdateTR = 0;
isSetUpdateVis = 0;
isSetCovReduce = 0;
isSetColor = 0;
isSetColorLambda = 0;

for i=1:2:numel(varargin)
    if ~isSetR && strcmpi(varargin{i},'r')
        
        R = varargin{i+1};
        R = R(:);
        isSetR = 1;
        
    elseif ~isSetT && strcmpi(varargin{i},'t')
        
        t = varargin{i+1};
        t = t(:);
        isSetT = 1;
        
    elseif ~isSetQ && strcmpi(varargin{i},'s')
        
        if isscalar(varargin{i+1})
            
            Q = repmat(varargin{i+1},K,1);
        else
            
            Q = varargin{i+1};
        end
                
        isSetQ = 1;
        Q = 1./Q;
        Q3 = repmat(Q, 1, 3);
    elseif ~isSetMaxNumIter && strcmpi(varargin{i},'maxnumiter')
        
        maxNumIter = varargin{i+1};
                
        isSetMaxNumIter = 1;
        
    elseif ~isSetEpsilon && strcmpi(varargin{i},'epsilon')
        
        epsilon = varargin{i+1};
                
        isSetEpsilon = 1;
        
    elseif ~isSetEpsilonST && strcmpi(varargin{i},'epsilonST')
        
        epsilonST = varargin{i+1};
                
        isSetEpsilonST = 1;    
        
    elseif ~isSetUpdatePriors && strcmpi(varargin{i},'updatepriors')
        
        updatePriors = varargin{i+1}; % don use the flag,it will affect subsequent parse
                
        isSetUpdatePriors = 1;
        
    elseif ~isSetUpdateTR && strcmpi(varargin{i}, 'updateTR')
        
        updateTR = varargin{i+1};
        
        isSetUpdateTR = 1;
       
    elseif ~isSetUpdateVis && strcmpi(varargin{i}, 'updateVis')
        
        updateVis = varargin{i+1};
        
        isSetUpdateVis = 1;
        
    elseif ~isSetCovReduce && strcmpi(varargin{i}, 'CovReduce')
        
        CovReduce = varargin{i+1};
        
        isSetCovReduce = 1;
        
    elseif ~isSetColor && strcmpi(varargin{i}, 'Color')
        
        Color = varargin{i+1};
        
        isSetColor = 1;
    elseif ~isSetColorLambda && strcmpi(varargin{i}, 'colorLambda')
        
        colorLambda = varargin{i+1};
        
        isSetColorLambda = 1;
        
        
        
    elseif ~isSetInitialPriorsOrGamma && strcmpi(varargin{i},'gamma')
        
        gamma = varargin{i+1};
                
        pk = repmat(1/(K*(gamma+1)),K,1);
        
        isSetInitialPriorsOrGamma = 1;
        
    elseif ~isSetInitialPriorsOrGamma && strcmpi(varargin{i},'initialpriors')
        
        if isscalar(varargin{i+1})
            
            pk = repmat(varargin{i+1},K,1);
        else
            
            pk = varargin{i+1};
            
        end
        
        gamma = (1-sum(pk))/sum(pk);
        
        isSetInitialPriorsOrGamma = 1;
        
    elseif strcmpi(varargin{i},'vis')    
        vis = varargin{i+1};
        isSetVis = 1;
    else
        
        if isSetInitialPriorsOrGamma
            
            error('Only one of the parameters ''initialPriors'' and ''gamma'' must be set.');
        else
            error('uknown option %s, or already set.',varargin{i});
        end
    end
end


% ======================================================================== %
%                   I N I T I A L I Z E   D E F A U L T S                  %
% ======================================================================== %

if ~isSetVis
    
    vis = {};
    for i=1:length(V)
       vis{end+1}=ones(1, K)./length(V);
    end
    
end

if ~isSetR
    
    R = repmat({eye(3)},M,1);
    
end

if ~isSetT
    
    t = cellfun(@(V) (-mean(V,2)+mean(X,2)),V,'uniformoutput',false);
    
end

% transformed sets based on inititial R & t (\phi(v) in the paper)
TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V,R,t,'uniformoutput',false);

[minXyZ,maxXyZ] = cellfun(@(x) deal(min(x,[],2),max(x,[],2)),[TV;X],'uniformoutput',false);
minXyZ = min(cat(2,minXyZ{:}),[],2);    
maxXyZ = max(cat(2,maxXyZ{:}),[],2);
QInit = repmat(1./(sqe(minXyZ,maxXyZ)),K,1);

if ~isSetQ
    
    [minXyZ,maxXyZ] = cellfun(@(x) deal(min(x,[],2),max(x,[],2)),[TV;X],'uniformoutput',false);
    
    minXyZ = min(cat(2,minXyZ{:}),[],2);
    
    maxXyZ = max(cat(2,maxXyZ{:}),[],2);
    
    Q = repmat(1./(sqe(minXyZ,maxXyZ)),K,1);
    Q3 = repmat(1./(sqe(minXyZ,maxXyZ)),K,3);
    Qcov = repmat({eye(3).*sqe(minXyZ,maxXyZ)}, K, 1);
end

if ~isSetMaxNumIter
    
    maxNumIter = 100;
    
end

if ~isSetEpsilon
    
    epsilon = 1e-6;
    
end

if ~isSetEpsilonST
    
    epsilonST = 0.05;
    
end

if ~isSetUpdatePriors
    
    updatePriors = 0;
    
end

if ~isSetUpdateTR
    
    updateTR = 1;
    
end

if ~isSetUpdateVis
    
    updateVis = 1;
    
end

if ~isSetCovReduce
    
    CovReduce = 1;
    
end

if ~isSetColor
    
    Color = {};
    
end

if ~isSetColorLambda
    
    colorLambda = 0;
    
end



if ~isSetInitialPriorsOrGamma
    
    gamma = 1/K;
    
    pk = repmat(1/(K+1),K,1);
    
end

% ======================================================================== %
%                                   E  M                                   %
% ======================================================================== %

% if requested, allocate an empty cell in maxNumIter dimension
if nargout > 6
    T = cell(M,2,maxNumIter);
    TAssign = cell(M,maxNumIter);
    TXQ = cell(2,maxNumIter);
end

% parameter h in the paper (this should be proportional to the volume that
% encompasses all the point sets). Above, we initially translate the sets around
% (0,0,0), and we compute accordingly the initial variances (and precisions)
%. Thus, we compute h in a similar way.
h = 2/mean(QInit); 

beta = gamma/(h*(gamma+1));
% beta = 4.1924e-04;
%keyboard

pk = pk'; % used as a row

% colorLambda=0.1;

for iter = 1:maxNumIter
    iter
    % POSTERIORS
    
    X_color = zeros(size(X));
    X_cnt = zeros(size(X, 2), 1);
    kdtree = vl_kdtreebuild(X);
    for ii=1:length(V)
        [index_single, ~] = vl_kdtreequery(kdtree, X, TV{ii}, 'NumNeighbors', 1, 'MaxComparisons', 10);
        colors = Color{ii}';
        sRMean = accumarray(index_single', colors(:,1), [], @mean);
        sGMean = accumarray(index_single', colors(:,2), [], @mean);
        sBMean = accumarray(index_single', colors(:,3), [], @mean);
    
        color_sub = [sRMean(:), sGMean(:), sBMean(:)]';
        X_color(:, 1:size(color_sub, 2)) = X_color(:, 1:size(color_sub, 2)) + color_sub;
        X_cnt(1:size(color_sub, 2)) = X_cnt(1:size(color_sub, 2)) + 1;
    end
    % sqe (squared differences between TV & X)
%     a = cellfun(@(TV, VC) sqe(TV,X),TV, Color, 'uniformoutput',false);
    
    %Color Version
    a = cellfun(@(TV, VC) (1-colorLambda)*sqe(TV,X)+colorLambda*(sqe(VC,X_color)),TV, Color, 'uniformoutput',false);
    
    %a3 = cellfun(@(TV) sqe2(TV,X),TV,'uniformoutput',false);
    
    % pk*S^-1.5*exp(-.5/S^2*||.||)
    a = cellfun(@(a) bsxfun(@times,pk.*(Q'.^1.5),exp(bsxfun(@times,-.5*Q',a))),a,'uniformoutput',false);
    %a3 = cellfun(@(a) bsxfun(@times,pk.*(Q'.^1.5),exp(bsxfun(@times,-.5*Q',a))),a,'uniformoutput',false);
    
%     a31 = cellfun(@(a) bsxfun(@times,pk.*(Q3(:,1)'.^1.5),exp(bsxfun(@times,-.5*Q3(:,1)',a))),a,'uniformoutput',false);
%     a32 = cellfun(@(a) bsxfun(@times,pk.*(Q3(:,2)'.^1.5),exp(bsxfun(@times,-.5*Q3(:,2)',a))),a,'uniformoutput',false);
%     a33 = cellfun(@(a) bsxfun(@times,pk.*(Q3(:,3)'.^1.5),exp(bsxfun(@times,-.5*Q3(:,3)',a))),a,'uniformoutput',false);
%     
%     a3 = cellfun(@(a1, a2, a3, vis) bsxfun(@times, ((a1+a2+a3) ./ 3), vis) , a31, a32, a33, vis', 'uniformoutput',false);

    % normalize
    a = cellfun(@(a, vis) bsxfun(@rdivide,a,sum(bsxfun(@times, a, vis),2)+beta),a, vis','uniformoutput',false);    
   
    
    %Update visibility term
%      if iter == 1
    if updateVis
%         epsilonST = 0.05;
        ajk = cellfun(@(a) sum(a), a, 'uniformoutput', false);
        ajkMat = cell2mat(ajk);

        [ st ] = genUniformDist( ajkMat, epsilonST );
        vis = genVisFromPeriod(st, length(vis), K, epsilonST);
    end
%      end
    
    % ------  weighted UMEYAMA ------ 
    
    lambda = cellfun(@(a) sum(a)',a,'uniformoutput',false); % 1 x K rows
%     lambda2 = cellfun(@(a) mean(a)',a,'uniformoutput',false); % 1 x K rows
    
    W = cellfun(@(V,a) bsxfun(@times,V*a,Q'),V,a,'uniformoutput',false);
    
    % weights, b
    b = cellfun(@(lambda) lambda.*Q,lambda,'uniformoutput',false);
    
    % mean of W
    mW = cellfun(@(W) sum(W,2),W,'uniformoutput',false);
    
    % mean of X
    mX = cellfun(@(b) X*b,b,'uniformoutput',false);
    
    % sumOfWeights
    sumOfWeights = cellfun(@ (lambda) dot(lambda,Q),lambda,'uniformoutput',false);
    
    % P
    P = cellfun(@(W,sumOfWeights,mW,mX) X*W' - mX*mW'/sumOfWeights, W,sumOfWeights,mW,mX,'uniformoutput',false);
    
    
    % SVD
    [uu,~,vv] = cellfun(@svd,P,'uniformoutput',false);

    
    if updateTR
    
        % compute R and check reflection
        R = cellfun(@(uu,vv) uu*diag([1 1 det(uu*vv)])*vv',uu,vv,'uniformoutput',false);

        % solve for t
        t = cellfun(@(mW,mX,R,sumOfWeights) (mX-R*mW)/sumOfWeights,mW,mX,R,sumOfWeights,'uniformoutput',false);
    
    end
    
    
    
    
    % transformed sets
    TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V,R,t,'uniformoutput',false);
    
    
    % UPDATE X

    den = sum(cell2mat(lambda'),2)'; % den is used for S's update as well
%     den2 = sum(cell2mat(lambda2'),2)'; % den is used for S's update as well

    X = cellfun(@(TV,a) TV*a,TV,a,'uniformoutput',false);
%     X = cellfun(@(TV,a) TV*a./size(a,1),TV,a,'uniformoutput',false);
%     X = cellfun(@(X) mean(X)',X,'uniformoutput',false); % 1 x K rows
    
    X = sum(cat(3,X{:}),3);
%     X = sum(cell2mat(X'),2)';
    
    X = bsxfun(@rdivide,X,den);
%     X = bsxfun(@rdivide,X,den2);
    
    % UPDATE S
    
    % denominators for each j
%      wnormes = cellfun(@(TV,a) sum(a.*sqe(TV,X)), TV,a, 'uniformoutput',false);
     
     %color version
     wnormes = cellfun(@(TV,a,VC) sum(a.*((1-colorLambda)*sqe(TV,X)+colorLambda*(sqe(VC,X_color)))), TV,a,Color, 'uniformoutput',false);
     
     
%     wnormes2 = cellfun(@(TV,a) mean(a.*sqe(TV,X)), TV,a, 'uniformoutput',false);
%     wnormes2 = sum(cell2mat(wnormes2'),2)';

% tic;
%     wnormes31 = cellfun(@(TV,a) mean(a.*sqe(TV(1,:),X(1,:))), TV,a, 'uniformoutput',false);
%     wnormes32 = cellfun(@(TV,a) mean(a.*sqe(TV(2,:),X(2,:))), TV,a, 'uniformoutput',false);
%     wnormes33 = cellfun(@(TV,a) mean(a.*sqe(TV(3,:),X(3,:))), TV,a, 'uniformoutput',false);
% toc;

     Q = transpose(3*den ./ (sum(cat(3,wnormes{:}),3) + 3*den*epsilon));
%     Q = transpose(3*den2 ./ (sum(cat(3,wnormes2{:}),3) + 3*den2*epsilon));
    
%     Q3 = [transpose(3*den2 ./ (sum(cat(3,wnormes31{:}),3) + 3*den2*epsilon)) ...
%              transpose(3*den2 ./ (sum(cat(3,wnormes32{:}),3) + 3*den2*epsilon)) ...
%                 transpose(3*den2 ./ (sum(cat(3,wnormes33{:}),3) + 3*den2*epsilon))];
    
    Q3 = repmat(Q, 1, 3);
    
    % UPDATE pk
    
    if updatePriors
        
%         pk = den / ((gamma+1)*sum(den));
        pk = den2 / ((gamma+1)*sum(den2));
        
    end
    
    
    % populate T
    if nargout > 6
        T(:,1,iter) = R;
        
        T(:,2,iter) = t;
        
        maxIdx = cell(length(a),1);

        for ii=1:length(a)
            [m, i] = max(a{ii}');
            maxIdx{ii} = i;
        end
        
        TAssign(:, iter) = maxIdx;
        
        TXQ{1,iter} = X;
        TXQ{2,iter} = Q;
    end
    
end


% return variances
if nargout > 3
    
    Q = 1./Q;
    
end



