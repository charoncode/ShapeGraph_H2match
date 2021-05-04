% get_conn_comp.m: function that separates a shape graph into its component curves, i.e., a disjoint union of open curves or 
% a single closed curve.
%
% Input:
% shape        : structure containing shape graph
%
% Output:
% conn_comp    : cell array of structures for each component curve of the shape graph
% updated_shape: structure containing shape graph updated with new ordering of vertices, edges and edge weights
% cut_vertices : structure containing cut vertices and their degrees


function [conn_comp, updated_shape, cut_vertices] = get_conn_comp(shape)

% Set default signal, edge weights and topology for the shape graph
if ~isfield(shape,'f')
    shape.f = zeros(length(shape.x),1);
end

if ~isfield(shape,'rho')
    shape.rho = zeros(length(shape.G),1);
end

if ~isfield(shape,'topology')
    shape.topology = 'general';
end

% Check if shape graph is a single closed curve, otherwise apply Tarjan's algorithm
if size(shape.x,1) == size(shape.G,1)
    
    conn_comp = cell(1);
    conn_comp{1}.x = shape.x;
    conn_comp{1}.G = shape.G;
    conn_comp{1}.f = shape.f;
    conn_comp{1}.rho = shape.rho;
    conn_comp{1}.topology = 'closed';
    
    updated_shape = shape;   
    cut_vertices.x = [];
    cut_vertices.deg = [];
    
else
 
    % Construct weighted directed graph object
    g = digraph(shape.G(:,1), shape.G(:,2), shape.rho);

    % Update edge list and edge weights after reordering of vertices of the shape graph
    shape.G = g.Edges.EndNodes;
    shape.rho = g.Edges.Weight;

    % Calculate degree = outdegree + indegree of each vertex
    outdeg = outdegree(g); 
    indeg = indegree(g);
    deg = outdeg + indeg;

    % Declare vertices with degree >=3 as cut vertices
    cut_ind = find(deg >= 3);
    cut_vertices.x = shape.x(cut_ind,:);  % coordinates of cut vertices
    cut_vertices.deg = deg(cut_ind);  % degree of cut vertices

    % Replicate cut vertices according to their degrees and update connectivity matrix of shape graph
    if ~isempty(cut_ind)

         for c = 1:length(cut_ind)  % for each cut vertex

            % find adjacent edges to the cut vertex   
            outgoing_edges = find(shape.G(:,1) == cut_ind(c));  % outgoing edges    
            incoming_edges = find(shape.G(:,2) == cut_ind(c));  % incoming edges       
            adj_edges = sort([outgoing_edges ; incoming_edges]);

            % update all but one of these adjacent edges as follows...
            updated_edges = adj_edges(1:cut_vertices.deg(c)-1);  % indices of edges that will be updated

            for l = 1:length(updated_edges)

                % replicate the cut vertex and update signal
                shape.x = [shape.x ; cut_vertices.x(c,:) ]; 
                shape.f = [shape.f ; 0]; 

                % determine if the updated edge is outgoing or incoming
                outgo = find(shape.G(updated_edges(l),:) == cut_ind(c)); 

                % update edge by creating a connection to the replicated cut vertex instead of the original vertex
                if outgo == 1 
                    shape.G(updated_edges(l),1) = length(shape.x);  % if outgoing
                else 
                    shape.G(updated_edges(l),2) = length(shape.x);  % if incoming
                end

            end  

        end

    end

    % Apply Tarjan's algorithm to break shape graph into its components, i.e., into disjoint union of open curves 
    [conn_comp, ~] = tarjan(shape);

    % Update new shape struct
    updated_shape.x = [];
    updated_shape.G = [];
    updated_shape.f = [];
    updated_shape.rho = [];
    updated_shape.topology = 'general';   
    num_vertices = 0;  % counting variable for number of vertices in new_shape

    for k = 1:length(conn_comp)

        updated_shape.x = [updated_shape.x ; conn_comp{k}.x];  % vertices
        updated_shape.G = [updated_shape.G ; num_vertices + conn_comp{k}.G];  % edges
        updated_shape.f = [updated_shape.f ; conn_comp{k}.f];  % signal
        updated_shape.rho = [updated_shape.rho ; conn_comp{k}.rho];  % edge_weights

        num_vertices = num_vertices + length(conn_comp{k}.x);

    end

end

end