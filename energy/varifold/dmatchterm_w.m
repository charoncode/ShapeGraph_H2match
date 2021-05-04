function [dxg, drhog] =  dmatchterm_w(fshape1,fshape2,objfun)
% g = dmatchterm_w(fshape1,fshape2,objfun) computes the derivative of the distance 
% function between fshape1 and fshape2 computed given by various currents 
% or varifold norms. Options are given by the structure objfun.

% Input:
%   fshape1 : a fshape structure
%   fshape2 : a fshape structure
%   objfun : structure containing the parameters : ('distance'=='cur','var' or 'varexpo'),('kernel_size_geom',geometric kernel bandwidth), ('sigmaf',functional kernel bandwidth), (method == 'cuda' or 'matlab') 

%Output:
% dxg  : derivative wrt x (n x d array) 
% drhog: derivative wrt weights rho (n x 1 array) (or probably M x 1 array, where M is number of edges in fshape2)

% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

% Some checks
if (size(fshape1.x,2) ~= size(fshape2.x,2)) || (size(fshape1.G,2) ~= size(fshape2.G,2))
    error('fshapes should be in the same space')
end


switch objfun.distance
    case 'kernel'
        switch objfun.kernel_distance.distance
            case 'var' % unoriented varifold with binet kernel
                DG = @dfvarifoldnorm_binet_w;
                
            case 'cur' % current distance (oriented)
                DG = @dfcurrentnorm_w;
                
            otherwise
                DG = @dfshape_kernel_distance_w;
        end
        
%     case 'wasserstein'
%         DG = @dfshape_wasserstein_distance;
%         
%     case 'pt2pt' % l2 distance
%         DG = @dfpt2pt;
end


[dxg,~,drhog] = DG(fshape1,fshape2,objfun);

end
