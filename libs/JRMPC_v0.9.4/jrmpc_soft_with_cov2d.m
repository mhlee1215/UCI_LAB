function [R,t,X,Q3,a,pk,T, TAssign, TXQ, vis, XC, XN] = jrmpc_soft_with_cov2d(V,X,varargin)
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

% sqe = @(Y,X) pdist2(X', Y')'.^2;% sum(bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1])).^2,3);
% sqe2 = @(Y,X) bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1]));

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
isSetNormal = 0;
isSetNormalLambda = 0;
isSetColor = 0;
isSetColorLambda = 0;

XC = [];
XN = [];

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
    elseif strcmpi(varargin{i},'normal')    
        normal = varargin{i+1};
        XN = zeros(3, K);
        isSetNormal = 1;
    elseif strcmpi(varargin{i},'normalLambda')    
        normalLambda = varargin{i+1};
        isSetNormalLambda = 1;
    elseif strcmpi(varargin{i},'color')    
        color = varargin{i+1};
        XC = zeros(3, K);
        isSetColor = 1;
    elseif strcmpi(varargin{i},'colorLambda')    
        colorLambda = varargin{i+1};
        isSetColorLambda = 1;    
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
       vis{end+1}=ones(1, K);%./length(V);
    end
    
end

if ~isSetR
    
    R = repmat({eye(3)},M,1);
    
end

if ~isSetT
    
%     t = cellfun(@(V) (-mean(V,2)+mean(X,2)),V,'uniformoutput',false);
    t = cellfun(@(V) (-mean(V,2)+mean(X,2))*0,V,'uniformoutput',false);
    
end

% transformed sets based on inititial R & t (\phi(v) in the paper)
TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V,R,t,'uniformoutput',false);

if isSetNormal
    TN = cellfun(@(N,R) R*N,normal,R,'uniformoutput',false);
end

[minXyZ,maxXyZ] = cellfun(@(x) deal(min(x,[],2),max(x,[],2)),[TV;X],'uniformoutput',false);
minXyZ = min(cat(2,minXyZ{:}),[],2);    
maxXyZ = max(cat(2,maxXyZ{:}),[],2);
QInit = repmat(1./(sqe(minXyZ,maxXyZ)),K,1);

if ~isSetQ
    
    [minXyZ,maxXyZ] = cellfun(@(x) deal(min(x,[],2),max(x,[],2)),[TV;X],'uniformoutput',false);
    
    minXyZ = min(cat(2,minXyZ{:}),[],2);
    
    maxXyZ = max(cat(2,maxXyZ{:}),[],2);
    
    Q = repmat(1./(sqe(minXyZ,maxXyZ)),K,1);
    
    QV = 1./(sqe(minXyZ,maxXyZ));
    
    Q3 = repmat(inv(eye(3).*(QV)), 1, 1, K);%repmat(1./(sqe(minXyZ,maxXyZ)),K,3);
    Q3i = ones(K, 1).*(QV.^-(3));
    
    Q3N = repmat(eye(3).*100, 1, 1, K);%repmat(1./(sqe(minXyZ,maxXyZ)),K,3);
    Q3Ni = ones(K, 1).*100;
    
    Q3C = repmat(eye(3).*100, 1, 1, K);%repmat(1./(sqe(minXyZ,maxXyZ)),K,3);
    Q3Ci = ones(K, 1).*100;
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

if ~isSetNormal
    normal = {};
    
end

if ~isSetNormalLambda
    normalLambda = 0;
end

if ~isSetColor
    color = {};
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

pk = gpuArray(pk);
X = gpuArray(X);
TV = cm2cg(TV);%cellfun(@(a) gpuArray(a), TV, 'uniformoutput',false);
Q = gpuArray(Q);
Q3 = gpuArray(Q3);
Q3i = gpuArray(Q3i);

gg = gpuDevice();

vRatio = 1/(1-normalLambda-colorLambda);
cRatio = colorLambda/(1-normalLambda-colorLambda);
nRatio = normalLambda/(1-normalLambda-colorLambda);

progressbar2(0, 0);
progressbar2('JRMPC registration...', 'Sub progress...');
for iter = 1:maxNumIter
%     iter
    fprintf('%d / %d (memory : %d)\n', iter, maxNumIter, gg.AvailableMemory);
    
    progressbar2([], 1/4);
    % POSTERIORS
    
    % sqe (squared differences between TV & X)
    %%%%%%%%%%%%%%%%%%%%%
%     %1D version %%%%%%
%     a = cellfun(@(TV) sqe(TV,X),TV,'uniformoutput',false);
%     
%     % merge with other features
%     a = cellfun(@(a) vRatio.*a,a,'uniformoutput',false);
%     
%     if isSetNormal
%         an = cellfun(@(TN) sqe(TN,XN),TN,'uniformoutput',false);
%         a = cellfun(@(a, an) a+nRatio.*an,a, an,'uniformoutput',false);
%     end
%     if isSetColor
%         ac = cellfun(@(C) sqe(C,XC),color,'uniformoutput',false);
%         a = cellfun(@(a, ac) a+cRatio.*ac,a, ac,'uniformoutput',false);
%     end
    %1D version end %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%

    
    
    % pk*S^-1.5*exp(-.5/S^2*||.||)

    %3D Ver.
    a = cellfun(@(TV) bsxfun(@times,pk.*(Q3i'.^-0.5),exp(-.5*getQSet(TV,X,Q3))),TV,'uniformoutput',false);
    
    a = cellfun(@(TV) bsxfun(@times,pk.*(Q3i'.^-0.5),exp((-.5*getQSet(TV,X,Q3)))),TV,'uniformoutput',false);
    
%     if isSetNormal
%         an = cellfun(@(TN) bsxfun(@times,pk.*(Q3Ni'.^-0.5),exp(-.5*getQSet(TN,XN,Q3N))),TN,'uniformoutput',false);
%         a = cellfun(@(a, an) a+nRatio.*an,a, an,'uniformoutput',false);
%     end
%     if isSetColor
%         ac = cellfun(@(C) bsxfun(@times,pk.*(Q3Ci'.^-0.5),exp(-.5*getQSet(C,XC,Q3C))),color,'uniformoutput',false);
%         a = cellfun(@(a, ac) a+cRatio.*ac,a, ac,'uniformoutput',false);
%     end


    a = cellfun(@(a, vis) bsxfun(@times, a, vis), a, vis','uniformoutput',false);    
    % normalize
    a = cellfun(@(a, vis) bsxfun(@rdivide,a,sum(bsxfun(@times, a, vis),2)+beta),a, vis','uniformoutput',false);    
  
    
    %Update visibility term
%      if iter == 1
    if updateVis %&& false
%         epsilonST = 0.05;
        ajk = cellfun(@(a) sum(a), a, 'uniformoutput', false);
        ajkMat = cell2mat(cg2cm(ajk));

        [ st ] = genUniformDist( ajkMat, epsilonST );
        vis = genVisFromPeriod(st, length(vis), K, epsilonST);
    end
%      end
    

    progressbar2([], 2/4);
    % ------  weighted UMEYAMA ------ 
    lambda = cellfun(@(a) sum(a)',a,'uniformoutput',false); % 1 x K rows
    
    W = cellfun(@(V,a) bsxfun(@times,V*a,(Q)'),V,a,'uniformoutput',false);
    
    % weights, b
    b = cellfun(@(lambda) lambda.*(Q),lambda,'uniformoutput',false);
    
    % mean of W
    mW = cellfun(@(W) sum(W,2),W,'uniformoutput',false);
    
    % mean of X
    mX = cellfun(@(b) X*b,b,'uniformoutput',false);
    
    % sumOfWeights
    sumOfWeights = cellfun(@ (lambda) dot(lambda,(Q)),lambda,'uniformoutput',false);
    
    % P
    P = cellfun(@(W,sumOfWeights,mW,mX) X*W' - mX*mW'/sumOfWeights, W,sumOfWeights,mW,mX,'uniformoutput',false);
    
    
    % SVD
    [uu,~,vv] = cellfun(@svd,P,'uniformoutput',false);

    
    if updateTR% && false
    
        % compute R and check reflection
        R2 = cellfun(@(uu,vv) uu*diag([1 1 det(uu*vv)])*vv',uu,vv,'uniformoutput',false);

        % solve for t
        t2 = cellfun(@(mW,mX,R,sumOfWeights) (mX-R*mW)/sumOfWeights,mW,mX,R,sumOfWeights,'uniformoutput',false);
        
        R2 = cg2cm(R2);
        t2 = cg2cm(t2);
        
        [R, t] = cellfun(@(R2, t2) fcn_supressConstraint(R2, t2),R2, t2,'uniformoutput',false);
    end
    
    
    
    
    % transformed sets
    TV = cellfun(@(V,R,t) bsxfun(@plus,R*V,t),V,R,t,'uniformoutput',false);
    
    if isSetNormal
        TN = cellfun(@(N,R) R*N, normal,R,'uniformoutput',false);
    end
    
%     if isSetColor
%         
%     end
    
    
    % UPDATE X
    progressbar2([], 3/4);
    lambda = cg2cm(lambda);
    den = sum(cell2mat(lambda'),2)'; % den is used for S's update as well
%     den2 = sum(cell2mat(lambda2'),2)'; % den is used for S's update as well

    X = cellfun(@(TV,a) TV*a,TV,a,'uniformoutput',false);  
    X = sum(cat(3,X{:}),3);
    X = bsxfun(@rdivide,X,den);

    if isSetNormal
        XN = cellfun(@(TN,a) TN*a,TN,a,'uniformoutput',false);
        XN = sum(cat(3,XN{:}),3);
        XN = bsxfun(@rdivide,XN,den);
    end
    
    if isSetColor
        XC = cellfun(@(color,a) color*a,color,a,'uniformoutput',false);
        XC = sum(cat(3,XC{:}),3);
        XC = bsxfun(@rdivide,XC,den);
    end


    % UPDATE S
    progressbar2([], 4/4);
    % denominators for each j
    dim = 3;
     wnormes = cellfun(@(TV,a) sum(a.*sqe(TV,X)), TV,a, 'uniformoutput',false); %M x 2 cell, 1xK each
     wnormes = cellfun(@(w) vRatio .* w, wnormes, 'uniformoutput', false); %M x 2 cell, 1xK each
     
%      wnormes3D = {};
%      for i=1:length(wnormes)
%          www = wnormes{i};
%          wnormes3D_sub = gpuArray(zeros(1, K, dim*dim));
%          wnormes3D_sub(:,:,1) = www;
%          wnormes3D_sub(:,:,5) = www;
%          wnormes3D_sub(:,:,9) = www;
%          wnormes3D{i} = wnormes3D_sub;
%      end
     
     wnormes3D = cellfun(@(TV,a) sum(repmat(a, 1, 1, dim*dim).*cell2mat(reshape(sqe(TV,X, true), 1, 1, dim*dim))), TV,a, 'uniformoutput',false); %M x 2 cell, 1xKx3x3 each
     wnormes3D = cellfun(@(w) vRatio .* w, wnormes3D, 'uniformoutput', false); %M x 2 cell, 1xK each
     if isSetNormal
%          wnormesN = cellfun(@(TN,a) sum(a.*sqe(TN,XN)), TN,a, 'uniformoutput',false);
         wnormes3DN = cellfun(@(TN,a) sum(repmat(a, 1, 1, dim*dim).*cell2mat(reshape(sqe(TV,X, true), 1, 1, dim*dim))), TN,a, 'uniformoutput',false); %M x 2 cell, 1xKx3x3 each
         wnormes3D = cellfun(@(w, wN) w + nRatio.*wN, wnormes3D, wnormes3DN, 'uniformoutput', false);
     end
%      
     if isSetColor
%          wnormes3DV = cellfun(@(color,a) sum(a.*sqe(color,XC)), color,a, 'uniformoutput',false);
         wnormes3DC = cellfun(@(TN,a) sum(repmat(a, 1, 1, dim*dim).*cell2mat(reshape(sqe(TV,X, true), 1, 1, dim*dim))), TN,a, 'uniformoutput',false); %M x 2 cell, 1xKx3x3 each
         wnormes3D = cellfun(@(w, wC) w + cRatio.*wC, wnormes3D, wnormes3DC, 'uniformoutput', false);
     end
     
%     wnormes2 = cellfun(@(TV,a) mean(a.*sqe(TV,X)), TV,a, 'uniformoutput',false);
%     wnormes2 = sum(cell2mat(wnormes2'),2)';

% tic;
%     wnormes31 = cellfun(@(TV,a) mean(a.*sqe(TV(1,:),X(1,:))), TV,a, 'uniformoutput',false);
%     wnormes32 = cellfun(@(TV,a) mean(a.*sqe(TV(2,:),X(2,:))), TV,a, 'uniformoutput',false);
%     wnormes33 = cellfun(@(TV,a) mean(a.*sqe(TV(3,:),X(3,:))), TV,a, 'uniformoutput',false);
% toc;

     Q = transpose(3*den ./ (sum(cat(3,wnormes{:}),3) + 3*den*epsilon));
     Q2 = transpose((sum(cat(3,wnormes{:}),3) + 3*den*epsilon)./(3*den));
     %For 3D version
     Q3 = gpuArray(zeros(dim, dim, K));
     cnt = 1;
     wnormes3D_sum = sum(cat(4,wnormes3D{:}),4);
     for ii=1:dim
         for jj=1:dim
            % Sigma^2
            Q3(ii,jj,:) = transpose((wnormes3D_sum(:,:,cnt) + 3*den*epsilon)./(3*den));
            cnt = cnt + 1;
         end
     end
     [Q3i] = q3inverse(Q3);
     
    
    % UPDATE pk
    if updatePriors        
        pk = den / ((gamma+1)*sum(den));
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
        TXQ{2,iter} = Q3;
    end
    
    progressbar2(iter / maxNumIter, []);
    
end


% return variances
if nargout > 3
    
%     Q = 1./Q;
    
end

%R,t,X,Q,a,pk,T
R = cg2cm(R);
t = cg2cm(t);
X = gather(X);
if ~isempty(XC)
    XC = gather(XC);
end
if ~isempty(XN)
    XN = gather(XN);
end
% Q = gather(Q);
a = cg2cm(a);
%pk;
T = cg2cm(T);

end

function [cm] = cg2cm(cg)
    cm = cellfun(@(a) gather(a), cg, 'uniformoutput',false);
    clear cg;
end

function [cg] = cm2cg(cm)
    cg = cellfun(@(a) gpuArray(a), cm, 'uniformoutput',false);
    clear cm;
end

function qSet = getQSet(TV, X, Q3, Q3i2)
    qSet = gpuArray(zeros(size(TV, 2), size(X, 2)));
    for i=1:size(X, 2)
        qSet(:,i) = getQ(TV, X(:,i), Q3(:,:,i));
    end
end

function q = getQ(X, mu, Sigma)
    X = bsxfun(@minus,X,mu);
    [U,p]= chol(Sigma);
    if p ~= 0
        Sigma
        error('ERROR: Sigma is not PD.');
        
    end
    Q = U'\X;
    q = dot(Q,Q,1);  % quadratic term (M distance)
end

function [q3i2] = q3inverse(q3)
    q3i2 = gpuArray(zeros(size(q3, 3), 1));
    for i=1:size(q3, 3)
        Sigma = q3(:,:,i);
        [U,p]= chol(Sigma);
        if p ~= 0
            error('ERROR: Sigma is not PD.');
        end
        q3i2(i) = prod(diag(U)).^2;
    end    
end

function [ee] = sqe(Y, X, is3D)% = @(Y,X) sum(bsxfun(@minus,permute(Y,[2 3 1]),permute(X,[3 2 1])).^2,3);

     YY = gather(Y); clear Y;
     XX = gather(X); clear X;
     
     if nargin > 2%exist('3D', 'var')
         dimension = size(XX, 1);
         ee = {};%zeros(size(YY, 2), size(XX, 2), dimension, dimension);
         ee2 = {};
         for i=1:dimension
%             ee2{i} = pdist2(XX(i,:)', YY(i,:)')';
            ee2{i} = (repmat(XX(i,:), size(YY(i,:),2), 1) - repmat(YY(i,:)', 1, size(XX(i, :), 2)));
         end
         for i=1:dimension
             for j=1:dimension
                  ee{i,j} = ee2{i}.*ee2{j};
             end
         end
     else
        ee = pdist2(XX', YY')'.^2;    
     end
end
