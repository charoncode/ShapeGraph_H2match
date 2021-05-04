% finishSplineSetup.m: function that adds knot sequences and collocation matrices for quadrature knots to splineData.
%
% Input: 
% splineData: structure containing essential B-spline parameters
%
% Output:
% splineData updated with the fields:
%   a) splineData.knotsS, splineData.innerKnotsS - knot sequences for spatial discretization
%   b) splineData.interpolS, splineData.noInterpolS - interpolation data for spatial discretization
%   c) splineData.knotsT, splineData.innerknotsT - knot sequences for time discretization
%   d) splineData.quadData, splineData.quadDataTensor - collocation matrices for quadrature knots for B-splines/tensor product B-splines


function splineData = finishSplineSetup(splineData)

% Construct knots
splineData = constructKnots(splineData);

% Construct collocation matrices
splineData = setupQuadData(splineData);
              
end


