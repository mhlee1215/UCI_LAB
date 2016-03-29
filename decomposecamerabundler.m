function [ K, R, t, C ] = decomposecamerabundler(P, k_diag)

if nargin == 1
    k_diag = [-1 -1 1];
end

M = (P(:,1:3));
C = -M\P(:,4);

[K, R] = rq(M);

% make diagonal of K = [- - +];
T = diag(sign(diag(K))./k_diag');
K = K * T;
R = T * R; % (T is its own inverse)

% Check negative det
if det(R) < 0
    warning('Determinant of rotation is negative.');
    R = -R;
end

% Get translation
t = -R*C;

% Normalize K (so last element is 1)
K = K/K(end);

end

