% edge_adjacency.m: function that computes the oriented edge adjacency matrix of an oriented graph.
%
% Input:
% G: edge list for a graph
%
% Output:
% E_adj: edge adjacency matrix for the graph

function [E_adj] = edge_adjacency(G)

% Preallocate sparse array for edge adjacency matrix
E_adj = spalloc(size(G,1), size(G,1), 3*size(G,1));

% Contruct edge adjacency matrix
for k = 1:size(G,1)  % for each edge of the graph
    
    % find destination vertex of that edge
    v2 = G(k,2);  
   
    % find indices of vertices connected to that destination vertex
    I1 = find(G(:,1) == v2);  
    I2 = find(G(:,2) == v2); 
    I2 = setdiff(I2, k);
    I = [I1 ; I2];
    
    % Update edge adjacency matrix
    for p = 1:length(I)
       E_adj(k, I(p)) = 1; 
    end
    
    % Symmetrize edge adjacency matrix
    E_adj = max(E_adj, E_adj');
    
end
