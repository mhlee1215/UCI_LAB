
addpath(genpath('libs'));


Vertices = bsxfun(@times, rand(3, 100), [1 3 1]');
X = mean(Vertices,2);
Q = cov(Vertices');
figure; scatter3(Vertices(1,:), Vertices(2,:), Vertices(3,:));
hold on;
h1 = plot_gaussian_ellipsoid(X, diag(diag(Q)), mean(mean(Q))*10 );
set(h1,'facealpha',0.2);
axis equal;

Q2 = covDimReduction(Q);
figure; scatter3(Vertices(1,:), Vertices(2,:), Vertices(3,:));
hold on;
h1 = plot_gaussian_ellipsoid(X, Q2, mean(mean(Q))*10 );
set(h1,'facealpha',0.2);
axis equal;

%  figure;
%  for i=length(T)
%      curX = TXQ{1,i};
%      curQ = TXQ{2,i};
%      for j=1:size(curX, 2)
%         h1 = plot_gaussian_ellipsoid(curX(:,j), [curQ(j) 0 0 ; 0 curQ(j) 0 ; 0 0 curQ(j)]./max(curQ), curQ(j)/max(curQ)*0.03);
%      end
%  end
%  
%  set(h1,'facealpha',0.6);
%  view(129,36); set(gca,'proj','perspective'); grid on; 
%  grid on; axis equal; axis tight;