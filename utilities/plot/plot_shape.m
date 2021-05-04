% plot_shape.m: function that plots a shape graph while respecting its graph stucture.
%
% Input: 
% shape: structure containing shape graph
% colorstyle: string specifying color/style of plot (optional)
% linewidth: specifies linewidth (optional)
%
% Output: plot of shape graph
%
% Note: This function handles both curves in R^2 and R^3.

function plot_shape(shape,varargin)

% Dimension
[~,d] = size(shape.x);

% Set colorstyle to blue and linewidth to 2 by default
colorstyle = 'b-';
linewidth = 2;

% Use user-specified parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii}, 'char'))
        switch (lower(varargin{ii}))
            case 'colorstyle'
                ii = ii + 1;
                colorstyle = varargin{ii};
            case 'linewidth'
                ii = ii + 1;
                linewidth = varargin{ii};
        end
    end
    ii = ii + 1;
end

% Plot shape
for k = 1:size(shape.G,1)
    
    hold on   
    if d == 2
        
        plot([shape.x(shape.G(k,1),1) shape.x(shape.G(k,2),1)],...
             [shape.x(shape.G(k,1),2) shape.x(shape.G(k,2),2)],...
             colorstyle,'LineWidth',linewidth) 
         
    elseif d == 3
        
        plot3([shape.x(shape.G(k,1),1) shape.x(shape.G(k,2),1)],...
              [shape.x(shape.G(k,1),2) shape.x(shape.G(k,2),2)],...
              [shape.x(shape.G(k,1),3) shape.x(shape.G(k,2),3)],...
              colorstyle,'LineWidth',linewidth)
        grid on
        view(d)
          
    end
    
end

end