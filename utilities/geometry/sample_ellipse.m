% sample_ellipse.m: function that uniformly samples N points from an ellipse with
%                   semi-major axis p and semi-minor axis q
%
% Input:
%   p: semi-major axis length
%   q: semi-minor axis length
%   N: number of points
%
% Output:
%   ellipse: vertices of the ellipse [Nx2 array]

function ellipse = sample_ellipse(p,q,N)

ellipse = zeros(N,2);

for i = 1:N
    
    ellipse(i,:) = [p*cos(((i-1)*2*pi)/N), q*sin(((i-1)*2*pi)/N)];
    
end

end


