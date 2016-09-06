function [ TV ] = drawTransformation( V, T, params )
%DRAWTRANSFORMATION Summary of this function goes here
%   Detailed explanation goes here

% marker, markerSize, clrmap
% h=figure;
if ~isfield(params, 'type')
    type = 1;
else
    type = params.type;
end

if type == 2
    if ~isfield(params, 'GMM_color')
        GMM_color = bsxfun(@plus, rand(params.K, 3)./2, [0.3 0.3 0.3]);
    else
        GMM_color = params.GMM_color;
    end
end

interval = params.interval;



M = length(V);
idx = transpose(1:M);
strIdx = arrayfun(@(j) sprintf('view %d',j),idx,'uniformoutput',false);
clrmap = {[1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0] ...
    ; [1 1 0]; [1 0 1] ; [0 1 1] ; [1 0 0] ; [0 1 0] ; [0 0 1] ; [1 1 1] ; [0 0 0] ; [0.3 0.3 0.3] ; [.5 0.9 0.9] ; [0.8 1 1]; [0 0 0]};
clrmap = clrmap(1:M);

markerSize = {};
for ii=1:M
    markerSize{end+1} = 5+ii*2;
end
markerSize = markerSize';
markerSize = markerSize(1:M);

marker_set = {'s', 'x', '.', '^'};
marker = {};
for i=1:M
    marker{end+1} = marker_set{mod(i, 4)+1};
end
marker = marker';



TXQ = params.TXQ;

maxNumIter = size(T, 3);
    for iter = 0:interval:maxNumIter-1
        % apply transformation of iteration : iter
        if iter == 0
            TV = V;
        else
            TV = cellfun(@(V,R_iter,t_iter) bsxfun(@plus,R_iter*V,t_iter),V,T(:,1,iter),T(:,2,iter),'uniformoutput',false);
        end

        clf;

        hold on, grid on

        title(sprintf('Registration of the sets after %d iteration(s).\n',iter),'fontweight','bold','fontsize',12);

        if type == 1
%         hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize.*1.5,[0 0 0],'filled'), TV, clrmap, marker, markerSize, 'UniformOutput', false);
            hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,clrmap,'filled'), TV, clrmap, marker, markerSize, 'UniformOutput', false);
            scatter3(TXQ{1, iter+1}(1,:)', TXQ{1, iter+1}(2,:)', TXQ{1, iter+1}(3,:)', 8, 'k', 'filled');
%         hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,clrmap,marker), TV, clrmap, marker, markerSize, 'UniformOutput', false);
%         legend(strIdx{:});
        elseif type == 2
            Assigned = params.Assigned(:,iter+1);
            hg2 = cellfun(@(TV,Assigned,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,GMM_color(Assigned,:),'filled'), TV, Assigned, markerSize(1:length(TV)), 'UniformOutput', false);
            scatter3(TXQ{1, iter+1}(1,:)', TXQ{1, iter+1}(2,:)', TXQ{1, iter+1}(3,:)', 8, 'k', 'filled');
        end

        axis equal;
    %     set(gca, 'position',get(gca,'position')+[+580 0 0 0]);

        % iteration 1 locks the axes of subsequent plots
        if iter == 0
           XLim = get(gca,'XLim');

           YLim = get(gca,'YLim');

           Zlim = get(gca,'ZLim');

    %        set(gca,'fontweight','bold','children',hg2);
             set(gca,'fontweight','bold');
        else
    %        set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold','children',hg2); 
             set(gca,'XLim',XLim,'YLim',YLim,'ZLim',Zlim,'fontweight','bold'); 
        end

        if isfield(params, 'view') && ~isempty(params.view)
         view(params.view);
        end

        hold off
        
        if isfield(params, 'pause')
            pause(params.pause);
        else
            pause(.12);    
        end
        
    end


end

