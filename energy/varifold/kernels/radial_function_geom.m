function r=radial_function_geom(x,derivative_order,objfun)
% r=RADIAL_FUNCTION_GEOM(x,opt,sig) implements the kernel and its derivatives :
%
%
% Input :
%   x : a matrix
%   derivative order: integer 0 or 1 typically
%   objfun : a structure with a field 'kernel_geom' containing a vector with kernel bandwidths
%
% Output :
%   r : matrix of the size of x
%
% See also : radial_function_sphere, radial_function_signal
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2017)

r=zeros(size(x));



switch lower(objfun.kernel_geom)
    case 'gaussian'
        
        for l=objfun.kernel_size_geom
            if derivative_order==0
                r=r + exp(-x/l^2);
            elseif derivative_order==1
                r=r -exp(-x/l^2)/l^2;
            end
        end
        
    case 'cauchy'
        
        for l=objfun.kernel_size_geom
            if derivative_order==0
                r=r + 1 ./ (1 + (x/l^2));
            elseif derivative_order==1
                r=r -1 ./ (l^2 * (1 + (x/l^2)) .^2);
            end
        end

    case 'compact'

        for l=objfun.kernel_size_geom
           t = (x/l^2);
           if t < 1
             if derivative_order==0
                r=r+(1 - t ) .^2;
              elseif derivative_order==1
                r=r -2 * (1- t);
            end
          end
        end
        
    case 'energy_distance'
        
        if derivative_order==0
            r=-sqrt(x);
        elseif derivative_order==1
            r=zeros(size(x));
            ind=find(x>eps);
            r(ind)=-1./(2*sqrt(x(ind)));
        end
        
                
end

end
