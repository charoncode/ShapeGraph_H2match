function g = matchterm_w(fshape1,fshape2,objfun)

% g = matchterm_w(fshape1,fshape2,objfun) computes the distance between
% fshape1 and fshape2 wrt various currents or varifold norms, with weights 
% placed on fshape2. Options are given by the structure objfun.

% Input:
%   fshape1 : fshape structure
%   fshape2 : fshape structure
%   objfun  : structure containing the parameters : ('distance'=='cur','var' or 'varexpo'),('kernel_size_geom',geometric kernel bandwidth), ('sigmaf',functional kernel bandwidth), (method == 'cuda' or 'matlab') 

% Output:
%   g = distance between fshape1 and fshape2 (real number)
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

% Some checks
if (size(fshape1.x,2) ~= size(fshape2.x,2)) || (size(fshape1.G,2) ~= size(fshape2.G,2))
    error('fshapes should be in the same space')
end

switch objfun.distance
    case 'kernel'
        switch objfun.kernel_distance.distance
            case 'var' % unoriented varifold with binet kernel
                G = @fvarifoldnorm_binet_w;
                
            case 'cur' % current distance (oriented)
                G = @fcurrentnorm_w;
                
            otherwise
                G = @fshape_kernel_distance_w;
        end
        
end


g = G(fshape1,fshape2,objfun);

end

