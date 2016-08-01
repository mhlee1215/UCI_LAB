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
marker = params.marker;
markerSize = params.markerSize;
clrmap = params.clrmap;
strIdx = params.strIdx;

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
%         hg2 = cellfun(@(TV,clrmap,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,clrmap,marker), TV, clrmap, marker, markerSize, 'UniformOutput', false);
%         legend(strIdx{:});
        elseif type == 2
            Assigned = params.Assigned(:,iter+1);
            hg2 = cellfun(@(TV,Assigned,marker,markerSize) scatter3(TV(1,:),TV(2,:),TV(3,:),markerSize,GMM_color(Assigned,:),'filled'), TV, Assigned, marker, markerSize, 'UniformOutput', false);
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

