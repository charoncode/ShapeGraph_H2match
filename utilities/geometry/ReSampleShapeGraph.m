% ReSampleShapeGraph.m: function that resamples shape graphs.
%
% Input:
% shape: structure containing shape graph, with vertices and edge list stored in a struct
% n    : number of points for resampling shape graph
%
% Output:
% resampled_shape: structure containing resampled shape graph

function resampled_shape = ReSampleShapeGraph(shape, n)

% Calculate scaling proportionality constant
kappa = n/length(shape.x);

% Get connected components
conn_comp = get_conn_comp(shape);
num_vertices_ar = 0;
num_vertices_br = 0;

% Resample each component and update its edge list, signal and edge weights
for k = 1:length(conn_comp) 
    
    % Extract and resample vertices
    comp_k = conn_comp{k}.x;
    conn_comp{k}.x = ReSampleCurve(comp_k',floor(kappa*length(comp_k))+1)';
    conn_comp{k}.x = conn_comp{k}.x(1:end-1,:); % remove last repeated point
    
    % Update number of vertices after resampling
    num_vertices_ar = num_vertices_br + length(conn_comp{k}.x); 
    
    % Check that the resampled shape has n vertices, and resample last component accordingly
    if (k == length(conn_comp)) && (num_vertices_ar ~= n)
        
        conn_comp{k}.x = ReSampleCurve(comp_k', (n-num_vertices_br)+1)';
        conn_comp{k}.x = conn_comp{k}.x(1:end-1,:); % remove last repeated point
        
    end
    
    % Update connectivity matrix, signal and edge weights
    conn_comp{k}.G = [1:length(conn_comp{k}.x)-1 ; 2:length(conn_comp{k}.x)]';
    conn_comp{k}.f = zeros(length(conn_comp{k}.x),1);
    conn_comp{k}.rho = zeros(length(conn_comp{k}.G),1);
    
    % Update number of vertices before next resampling
    num_vertices_br = num_vertices_br + length(conn_comp{k}.x);
    
end

% Update struct for resampled shape
resampled_shape.x = [];
resampled_shape.G = [];
resampled_shape.f = [];
resampled_shape.rho = [];
resampled_shape.topology = 'general';   
num_vertices = 0; % counting variable for number of vertices in resampled shape
   
for k = 1:length(conn_comp)
        
    resampled_shape.x = [resampled_shape.x ; conn_comp{k}.x]; % vertices
    resampled_shape.G = [resampled_shape.G ; num_vertices + conn_comp{k}.G]; % edges
    resampled_shape.f = [resampled_shape.f ; conn_comp{k}.f]; % signal
    resampled_shape.rho = [resampled_shape.rho ; conn_comp{k}.rho]; % edge weights
        
    num_vertices = num_vertices + length(conn_comp{k}.x);
    
end
    

end