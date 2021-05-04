% get_grad_avg.m: function that returns the gradient averaging matrix for a shape graph.
%
% Input:
%   shape             : structure containing shape graph
%   cut_vertices      : cut vertices of shape graph
%
% Output:
%   grad_avg_vertices : gradient averaging matrix wrt to the vertices
%   grad_avg_cPts     : gradient averaging matrix wrt to the control points


function [grad_avg_vertices, grad_avg_cPts] = get_grad_avg(shape, cut_vertices)

% Extract component curves
conn_comp = shape.connComp;

% Initialize gradient averaging matrices as identity matrices
grad_avg_vertices = speye(length(shape.x));
grad_avg_cPts = speye(length(shape.cPts));

if ~isempty(cut_vertices.x)
    
    for c = 1:size(cut_vertices.x,1)  % for each cut vertex
        
        % store cut vertex and its degree
        p = cut_vertices.x(c,:); 
        d_p = cut_vertices.deg(c); 
        
        % Update gradient averaging matrix wrt vertices and control points
        ind_vertices = [];  % vertex - cut vertex list
        ind_cPts = [];  % control points - cut vertex list
        l_vertices = 0;  % counting variable for vertex list
        l_cPts = 0;  % counting variable for control point list
        
        for k = 1:length(conn_comp)  % for each component curve
            
            % find if and where cut vertex appears as a vertex of shape
            [~, q_vertices] = ismember(p, conn_comp{k}.x, 'rows'); 
            
            if q_vertices > 0  % if cut vertex appears
                
                ind_vertices = [ind_vertices ; l_vertices + q_vertices];  % update vertex index list
                
                % find corresponding index in list of control points
                if q_vertices == 1  % if first vertex in component curve
                    
                    ind_cPts = [ind_cPts ; l_cPts + 1];  % update control point index list
                    
                else  % if last vertex in component
                    ind_cPts = [ind_cPts ; l_cPts + length(conn_comp{k}.cPts)];
                        
                end
                
            end
            
            l_vertices = l_vertices + length(conn_comp{k}.x);  % update counting variable for vertex list
            l_cPts = l_cPts + length(conn_comp{k}.cPts);  % update counting variable for control point list
            
        end
        
        % Update averaging matrix
        for s = 1:length(ind_vertices)     
            for t = 1:length(ind_vertices)
                
                grad_avg_vertices(ind_vertices(s), ind_vertices(t)) = 1/d_p;
                grad_avg_cPts(ind_cPts(s), ind_cPts(t)) = 1/d_p;
                
            end  
        end  
        
    end
                   
end

end