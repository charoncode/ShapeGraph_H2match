% quartic_pen.m: function that computes the double-well penalty of a weight function (rho) defined on the edges of a 
% shape graph (c), and its gradient wrt the edge weights and the vertices of the shape graph.
%
% Input:
%   shape        : structure containing the shape graph
%   m1, m2       : desired minimizers of the double-well penalty (optional, default is m1=-1, m2=0 for the {0,1}-penalty)
%   pen_function : functional form of double-well penalty (optional, default is 'polynomial') [string]
%                  Options are:
%                  (i)  'polynomial' - smooth double-well polynomial penalty computed using Hermite interpolation
%                  (ii) 'hat' - non-smooth double-well penalty where valleys are separated by a 'hat'
%   normalization: normalization of the double-well penalty by edge lengths of shape (optional, default is 'none') [string]
%                  Options are:
%                  (i)  'none' - no normalization by corresponding edge lengths (all edge lengths assumed to be 1)   
%                  (ii) 'length' - normalize by corresponding edge lengths
%   h            : height of local maximum of double-well penalty (optional, default is h = 0.5) [scalar]
%
%
% Output:
%   pen          : value of double-well penalty evaluated entrywise at rho [Ex1 array]
%   drhopen      : gradient of double-well penalty wrt to rho [Ex1 array]
%   dcpen        : differential of normalized double-well penalty wrt to c [Ex2 array]

function [pen, drhopen, varargout] = quartic_pen(shape,varargin)

global objfunc

% Set default options
    if nargin < 6
        m1 = -1; 
        m2 = 0; 
        pen_function = 'polynomial'; 
        normalization = 'none'; 
        h = 0.5;
    
    else % otherwise use specified values
        m1 = varargin{1};
        m2 = varargin{2};
        pen_function = varargin{3};
        normalization = varargin{4};
        h = varargin{5};
    end
    
% Dimensions
E = length(shape.rho); % number of edges

% Set cut-off margins for clipping double-well penalty and desired gradient beyond these cutoff points
clip_pt = 0.25; clip_grad = 0.5;
pen = zeros(E,1);
drhopen = zeros(E,1);

% Compute double-well penalty and its gradient wrt edge weights 
switch lower(pen_function)
    
    % smooth polynomial double-well penalty - computed using Hermite interpolation
    case 'polynomial' 

    % Compute midpoint of desired minimizers (zeros of the penalty)
    m = (m1 + m2)/2;

    % Matrix of polynomial expression and its derivative evaluated at m1, m2, m
    X = [m1.^(0:5); m2.^(0:5); m.^(0:5);...
         0 (1:5).*(m1.^(0:4)); 0 (1:5).*(m2.^(0:4)); 0 (1:5).*(m.^(0:4))];

    % Value of polynomial double-well penalty (and its derivatives) at m1, m2, m
    y = [0; 0; h; 0; 0; 0];

    % Compute coefficients of polynomial penalty function
    C = X\y; 
    deg = length(C)-1;  % degree

    % Evaluate smooth polynomial double-well penalty and its derivative at rho
        for k = 1:E
    
            if shape.rho(k) < m1 - clip_pt  % left cutoff point
                pen(k) = ((m1 - clip_pt).^(0:deg))*C + -clip_grad*(shape.rho(k) - (m1 - clip_pt));
                drhopen(k) = -clip_grad; 
        
            elseif shape.rho(k) > m2 + clip_pt  % right cutoff point
                pen(k) = ((m2 + clip_pt).^(0:deg))*C + clip_grad*(shape.rho(k) - (m2 + clip_pt));
                drhopen(k) = clip_grad;        
        
            else
                pen(k) = (shape.rho(k).^(0:deg))*C;
                drhopen(k) = (1:5).*(shape.rho(k).^(0:4))*C(2:deg+1);
        
            end
    
        end
    
        
    % double-well penalty with hat
    case 'hat' 
        
        % Compute midpoint of desired minimizers (zeros of the penalty)
        m = (m1 + m2)/2;

        % Evaluate hat double-well penalty and its derivative wrt to edge weights
        for k = 1:E
    
            if shape.rho(k) < m1 - 4*clip_pt
                pen(k) = (1/2)*((-4*clip_pt)*(m1 - 4*clip_pt - m2))^2 - clip_grad*(shape.rho(k) - (m1 - 4*clip_pt));
                drhopen(k) =  -clip_grad;
        
            elseif shape.rho(k) > m2 + 4*clip_pt
                pen(k) = (1/2)*((m2 + 4*clip_pt - m1)*(4*clip_pt))^2 + clip_grad*(shape.rho(k) - (m2 + 4*clip_pt));
                drhopen(k) = clip_grad;
        
            elseif shape.rho(k) > m1 && shape.rho(k) <= m
                pen(k) = (h/(m - m1))*(shape.rho(k) - m1);
                drhopen(k) = h/(m - m1);

            elseif shape.rho(k) > m && shape.rho(k) < m2
                pen(k) = (h/(m - m2))*(shape.rho(k) - m2);
                drhopen(k) = h/(m - m2);   
         
            else
                pen(k) = (1/2)*((shape.rho(k) - m1).*(shape.rho(k) - m2)).^2;
                drhopen(k) = (shape.rho(k) - m1).*(shape.rho(k) - m2).*(2*shape.rho(k) - m1 - m2);
    
            end
    
        end
        
end

% Normalize double-well penalty by length of edges (if specified), and compute its gradient wrt to edge weights and vertices of the shape graph
switch lower(normalization)
    
    case 'length' 

        % Dimensions
        [N,d] = size(shape.x); dE = size(shape.G,2); m = factorial(dE-1);
                
        % Compute tangent vectors (edges) to shape graph
        switch lower(shape.topology)
            
            case {'closed','open'}
                c_prime = pVectors(shape.x,shape.G,'centered') / m; 
                
            otherwise 
                c_prime = pVectors(shape.x,shape.G,'forward') / m; 
                    
        end
        
        % Compute length of edges of shape graph
        L = vecnorm(c_prime,2,2); 
                
        % Normalize double-well penalty by length of edges
        pen = pen .* L;
        
        % Compute derivative of length normalized double-well penalty wrt edge weights
        drhopen = drhopen .* L;
                
        % Compute derivative of normalized double-well penalty wrt transformed source c
        switch lower(objfunc.edge_weight) 
            
            case {'both','source'} % perform only if optimizing over source edge weights
                
                dc_prime = repmat(pen,1,d) .* c_prime .* repmat(1./L,1,d);
                
                % Construct derivative of double well penalty wrt to c using dc_prime
                switch lower(shape.topology)
    
                    case 'open'
     
                        Conn = [1 1:(N-1) ; 2:N N]'; % assume ordering 1,..,N of the vertices and open curve
   
                        if d == 2   
                            dcpen = [accumarray(Conn(:),[-dc_prime(1,1) ; -dc_prime(2:(N-1),1)/2 ; -dc_prime(N,1) ; dc_prime(1,1) ; dc_prime(2:(N-1),1)/2 ; dc_prime(N,1)] ,[N,1],[],0),...
                                     accumarray(Conn(:),[-dc_prime(1,2) ; -dc_prime(2:(N-1),2)/2 ; -dc_prime(N,2) ; dc_prime(1,2) ; dc_prime(2:(N-1),2)/2 ; dc_prime(N,2)] ,[N,1],[],0) ];
                        
                        elseif d == 3
                            dcpen = [accumarray(Conn(:),[-dc_prime(1,1) ; -dc_prime(2:(N-1),1)/2 ; -dc_prime(N,1) ; dc_prime(1,1) ; dc_prime(2:(N-1),1)/2 ; dc_prime(N,1)] ,[N,1],[],0),...
                                     accumarray(Conn(:),[-dc_prime(1,2) ; -dc_prime(2:(N-1),2)/2 ; -dc_prime(N,2) ; dc_prime(1,2) ; dc_prime(2:(N-1),2)/2 ; dc_prime(N,2)] ,[N,1],[],0),...
                                     accumarray(Conn(:),[-dc_prime(1,3) ; -dc_prime(2:(N-1),3)/2 ; -dc_prime(N,3) ; dc_prime(1,3) ; dc_prime(2:(N-1),3)/2 ; dc_prime(N,3)] ,[N,1],[],0),];    
                        end   
                            
                    case 'closed'
 
                        Conn = [N 1:(N-1) ; 2:N 1]';
                        
                        if d == 2
                            dcpen = [accumarray(Conn(:),[-dc_prime(:,1)/2 ; dc_prime(:,1)/2] ,[N,1],[],0),...
                                     accumarray(Conn(:),[-dc_prime(:,2)/2 ; dc_prime(:,2)/2] ,[N,1],[],0) ];
                                
                        elseif d == 3
                            dcpen = [accumarray(Conn(:),[-dc_prime(:,1)/2 ; dc_prime(:,1)/2] ,[N,1],[],0),...
                                     accumarray(Conn(:),[-dc_prime(:,2)/2 ; dc_prime(:,2)/2] ,[N,1],[],0),...
                                     accumarray(Conn(:),[-dc_prime(:,3)/2 ; dc_prime(:,3)/2] ,[N,1],[],0)];
                        end
    
                    otherwise

                        if d == 2
                            dcpen = [accumarray(shape.G(:),[-dc_prime(:,1); dc_prime(:,1)] ,[N,1],[],0),...
                                     accumarray(shape.G(:),[-dc_prime(:,2); dc_prime(:,2)] ,[N,1],[],0) ];
                        
                        elseif d == 3
                            dcpen = [accumarray(shape.G(:),[-dc_prime(:,1); dc_prime(:,1)] ,[N,1],[],0),...
                                     accumarray(shape.G(:),[-dc_prime(:,2); dc_prime(:,2)] ,[N,1],[],0),...
                                     accumarray(shape.G(:),[-dc_prime(:,3); dc_prime(:,3)] ,[N,1],[],0)];
                        end

                end
            
        % Return output 
        varargout{1} = dcpen;
                
        end       
                
end
        
end

