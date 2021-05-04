% remove_segment.m: function which removes at most two disjoint segments from a shape graph.
%
% Input:
% shape    : structure containing shape graph
% ind_list : cell array of indices indicating which vertices are to be removed
% topology : string declaring whether curve is open or closed (default is closed)
%
% Output:
% new_shape: structure containing vertices and connectivity matrix of the shape graph after 
% removal of specified vertices
%
% Note: This function handles shape graphs in R^2 and R^3.

function new_shape = remove_segment(shape, ind_list, varargin)

% Set defaults
if nargin < 3
    varargin{1} = 'closed';
end

% Initial coordinates of new shape
new_shape.x = shape.x;

% Remove segments
for j = 1:length(ind_list)
    
    % For each segment to be removed...
    ind = ind_list{j};
    
    % Remove vertices and redefine edge list
    new_shape.x = new_shape.x(setdiff(1:size(new_shape.x,1), ind),:);
    
    if j == 1
        
        switch varargin{1}
            
            case 'closed'
                
                new_shape.G = [1:ind(1)-2 ind(1):size(new_shape.x,1) ; ...
                               2:ind(1)-1 ind(1)+1:size(new_shape.x,1) 1]';
                           
            case 'open'
                
                new_shape.G = [1:ind(1)-2 ind(1):size(new_shape.x,1)-1; ...
                               2:ind(1)-1 ind(1)+1:size(new_shape.x,1)]';
                           
        end
               
    else
        
        switch varargin{1}
            
            case 'closed'
                
                new_shape.G = [1:ind(1)-2 ind(1):size(new_shape.x,1)-1 ; ...
                               2:ind(1)-1 ind(1)+1:size(new_shape.x,1) ]';
                            
            case 'open'
                
                new_shape.G = [1:ind(1)-2 ind(1):size(new_shape.x,1)-2 ; ...
                               2:ind(1)-1 ind(1)+1:size(new_shape.x,1)-1 ]';
                           
        end
                   
    end
       
    % Reorder vertices and update edge list
    [order, H] = toposort(digraph(new_shape.G(:,1), new_shape.G(:,2)));
    new_shape.x = new_shape.x(order,:);
    new_shape.G = H.Edges.EndNodes;
    
end

end
