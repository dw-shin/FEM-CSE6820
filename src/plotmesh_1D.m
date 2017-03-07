function plotmesh_1D(c4n, n4e, ~)
%%
% plotmesh_1D    1 dimensional mesh plot
%    plotmesh_1D(c4n,n4e) displays the interval defined in the 2-by-N
%    matrix n4e as a mesh. Each column of n4e contains indices into the
%    left-end and the right-end points of the corresponding element.
%
%    plotmesh_1D(c4n,n4e,onoff) displays the interval defined in the 2-by-N
%    matrix n4e as a mesh and all nodes for the approximate solution. 
%    In this case, onoff is an arbitrary number.
%
%     - Input
%      c4n    coordinates for nodes.
%             c4n is a N dimensional vector and it contains all
%             coordinates for nodes of the approximate solution.
%      n4e    nodes for elements.
%             n4e is a 2-by-M matrix. Each column of n4e containes indices 
%             into the left-end and the right-end points of the
%             corresponding element.

%%
figure('Position', [200 800 600 100]);
set(gcf, 'Color', 'w')

allNodes = unique(n4e(:));
plot(c4n(allNodes), zeros(1,length(allNodes)), 'k+-', 'LineWidth', 2, 'MarkerSize',14)
if nargin>2
    hold on
    plot(c4n, 0, 'ro', 'markersize', 7)
    
    % node number
    for k = 1:size(c4n,2)
        text(c4n(1, k), 0, num2str(k, 'v_{%d}'), ...
            'VerticalAlignment', 'bottom');
    end
    
    % element number
    for k = 1:size(n4e,2)
        mid = (c4n(:,n4e(1, k)) + c4n(:,n4e(2, k)))/2;
        text(mid(1), 0, num2str(k,'T_{%d}'), 'Color', 'b', ...
            'VerticalAlignment', 'top');
    end
    hold off
end
end