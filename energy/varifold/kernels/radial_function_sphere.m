function res=radial_function_sphere(x,derivative_order,objfun)
% r = RADIAL_FUNCTION_SPHERE(x,derivative_order,objfun) computes of rho(x^2/2) (ie opt==0), rho'(x^2/2) (ie opt==1) and
% rho"(x^2/2)  (ie opt==2)
%
% Input :
%   x : a matrix
%   derivative order: integer 0 or 1 typically
%   objfun : a structure with a field 'kernel_grass' containing a vector with kernel bandwidths
%
% Output :
%   r : matrix of the size of x
%
% See also : radial_function_geom, radial_function_signal
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, G. Nardi, A. Trouve (2012-2016)

res=zeros(size(x));

switch lower(objfun.kernel_grass)
    case 'linear'
        
        if (derivative_order==0)
            res=x;
        elseif (derivative_order==1)
            res=ones(size(x));
        end
        
    case 'gaussian_oriented'
        
        for l = objfun.kernel_size_grass
            if derivative_order==0
                res = res + exp( (-2 + 2*x) / l^2);
            elseif derivative_order==1
                res = res + 2 * exp( (-2 +2*x) / l^2) / l^2;
            end
        end
        
        
    case 'binet'
        
        if (derivative_order==0)
            res=x.^2;
        elseif (derivative_order==1)
            res=2*x;
        end
        
    case 'gaussian_unoriented'
        
        for l=objfun.kernel_size_grass
            if (derivative_order==0)
                res = res + exp( (-2 + 2 * x.^2) / l^2);
            elseif (derivative_order==1)
                res = res + 4 * x .* exp( (-2 + 2 *x.^2) /l^2) / l^2;
            end
        end
        
    otherwise
        
        if (derivative_order==0)
            res=objfun.distance_kernel(x);
        elseif (derivative_order==1)
            res=objfun.distance_kernel_der(x);
        end
        
end
end

