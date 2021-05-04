% diff_operator.m: function that computes the difference operator for each component curve of a shape graph and concatenates
%                  them into a block diagonal matrix. This difference operator is used to compute the TV norm of
%                  a weight function defined on the edges of that shape graph.
%
% Input:
%  shape: structure containing shape
%
% Output:
%   D  : difference operator, stored as sparse array

function D = diff_operator(shape)

% Get connected components
if ~isfield(shape,'connComp')    
    shape.connComp = get_conn_comp(shape);    
end

% Preallocate arrays for difference operator
D_k = cell(length(shape.connComp), 1);
D = [];

% Construct difference operator for each component curve
for k = 1:length(shape.connComp)  
    
    % Get number of edges for each component
    M_k = size(shape.connComp{k}.G, 1);
    
    % Create square matrix of -1s on the diagonal and 1s on the super diagonal
    D0 = diag(-ones(M_k, 1)) + diag(ones(M_k-1, 1), 1);
    
    if strcmp(shape.connComp{k}.topology,'open')  % if component curve is open
    
        % Remove last row of D0, as difference operator is an M_k-1 * M_k matrix with -1s on the
        % diagonal and 1s on the super diagonal
        D_k{k} = D0(1:end-1,:);
        
    else  % if component curve is closed
        
        % Difference operator is an M_k * M_k matrix of -1s on the diagonal and 1s on the super
        % diagonal, with the last row having a 1 in first column and -1 on the last column
        D_k{k} = D0;
        D_k{k}(end,1) = 1;
        
    end
    
    % Create difference operator D by concatenating the D_k matrices into a big block diagonal matrix
    D = blkdiag(D, D_k{k});
    
end

% Save difference operator as sparse array
D = sparse(D);

end
