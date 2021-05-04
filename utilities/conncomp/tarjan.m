% tarjan.m: function that separates oriented shape graphs in R^d into disjoint component curves.
%
% Input:
% shape        : structure containing shape graph
%
% Output:
% conn_comp    : cell array of structures for each component curve of the shape graph
% updated_shape: structure containing shape graph updated with new ordering of vertices, edges and edge weights


function [conn_comp, updated_shape] = tarjan(shape)

% Set default signal, edge weights and topology for the shape graph
if ~isfield(shape,'f')
    shape.f = zeros(length(shape.x), 1);
end

if ~isfield(shape,'rho')
    shape.rho = zeros(length(shape.G), 1);
end

if ~isfield(shape,'topology')
    shape.topology = 'general';
end

% Construct directed weighted graph object for the shape graph
g = digraph(shape.G(:,1), shape.G(:,2), shape.rho);

% Reorder vertices, edges and edge weights
[order, G] = toposort(g);
updated_shape.x = shape.x(order,:);
updated_shape.G = G.Edges.EndNodes;
updated_shape.f = shape.f(order,:);
updated_shape.rho = G.Edges.Weight;
updated_shape.topology = shape.topology;

% Use Tarjan's algorithm to extract component curves
[num_comp, vertex_labels] = graphconncomp(adjacency(G), 'Weak', true);

% Convert component curves into fshape structure
conn_comp = cell(num_comp, 1);

for k = 1:num_comp  % for each component curve...
    
    % List vertices belonging to that component curve
    v = find(vertex_labels == k); 
    
    % Construct subgraph object for that component curve
    comp = subgraph(G, v);
    
    % Store each component curve as fshape structure
    conn_comp{k}.x = updated_shape.x(v,:);  % vertices
    conn_comp{k}.G = comp.Edges.EndNodes;  % edges (connectivity matrix)
    conn_comp{k}.f = updated_shape.f(v,:);  % signal
    conn_comp{k}.rho = comp.Edges.Weight;  % edge weights
    
    if length(v) == size(conn_comp{k}.G,1)   
        conn_comp{k}.topology = 'closed';  % closed curve
        conn_comp{k}.G = [1:length(v) ; 2:length(v) 1]';  % update edge list 
        
    else
        conn_comp{k}.topology = 'open'; % open curve
    end
    
end

end