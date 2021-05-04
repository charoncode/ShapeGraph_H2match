function r=radial_function_signal(x,derivative_order,objfun)
% r=RADIAL_FUNCTION_SIGNAL(x,derivative_order,objfun) implements the kernel and its derivatives :
%
%
% Input :
%   x : a matrix
%   derivative order: integer 0 or 1 typically
%   objfun : a structure with a field 'kernel_signal' containing a vector with kernel bandwidths
%
% Output :
%   r : matrix of the size of x
%
% See also : radial_function_geom, radial_function_sphere
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2017)

r=zeros(size(x));



switch lower(objfun.kernel_signal)
    case 'gaussian'
        
        for l=objfun.kernel_size_signal
            if derivative_order==0
                r=r + exp(-x/l^2);
            elseif derivative_order==1
                r=r -exp(-x/l^2)/l^2;
            end
        end
        
    case 'cauchy'
        
        for l=objfun.kernel_size_signal
            if derivative_order==0
                r=r + 1 ./(1 + (x/l^2));
            elseif derivative_order==1
                r=r -1 ./ (l^2 * (1 + (x/l^2)) .^2);
            end
        end
        
        
end

end
