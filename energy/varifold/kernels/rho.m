function r=rho(x,opt,sig)
% r=RHO(x,opt,sig) implements the Gaussian kernel and its derivatives :
%
%       rho(t)=exp(-t/lam^2);
%
% Computation of rho(x^2/2) (ie opt==0), rho'(x^2/2) (ie opt==1) and
% rho"(x^2/2)  (ie opt==2) 
%
% Input :
%   x : a matrix
%   opt : a integer
%   sig : vector with kernel bandwidths
%
% Output :
%   r : matrix of the size of x
%
% See also : dsOptim.c, Kernel3Deven.m, *.cu files for speed improvement via mex files.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

r=zeros(size(x));
for l=sig
    
    if opt==0
        r=r + exp(-x/l^2);
    elseif opt==1
        r=r -exp(-x/l^2)/l^2;
    elseif opt==2
        r=r +exp(-x/l^2)/l^4;
    end
end

end
