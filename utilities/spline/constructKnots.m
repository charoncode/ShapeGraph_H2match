%% constructKnots
%
% Function constructs the knot sequences with the given parameters from
% splineData. The knot sequences are uniform in space and time, periodic in
% space and with full multiplicity at the boundary in time.
%
% interpolS is set automatically using N and nS.
%
% Input
%   splineData
%       Contains information about spline degree and number of control
%       points.
%
% Output
%   splineData
%       Same as input with knot sequences.
%
% Note: This function is part of the H2 metrics library (https://github.com/h2metrics/h2metrics).


function [ splineData ] = constructKnots( splineData )

% Rename control point variables
N = splineData.N;
nS = splineData.nS;
Nt = splineData.Nt;
nT = splineData.nT;
curveClosed = splineData.curveClosed;

% Contruct knots for spatial discretization of curves in the path of the transformed source
if curveClosed % if the curve is closed

    % Normalize, domain of definition is [0,2*pi]
    splineData.knotsS = ((-nS):(N+nS))/N*2*pi; 
    splineData.innerKnotsS = splineData.knotsS(nS+1:end-nS);

    splineData.noInterpolS = splineData.N;
    innerKnotsS = splineData.innerKnotsS;

        if mod(splineData.nS, 2) == 0 % See deBoor (p.287) for reasons.
            splineData.interpolS = innerKnotsS(1:end-1)' + 0.5*diff(innerKnotsS)';
        else
            splineData.interpolS = innerKnotsS(1:end-1)';
        end

else % if curve is open

    % Normalize, domain of definition is [0,2*pi]
    splineData.knotsS = [zeros(1,nS), linspace(0,2*pi,N-nS+1), 2*pi*ones(1,nS)];
    splineData.innerKnotsS = splineData.knotsS(nS+1:end-nS);

    splineData.noInterpolS = splineData.N;
    splineData.interpolS = aveknt(splineData.knotsS, nS+1)';

end

% Construct knots for time discretization
splineData.knotsT = [ zeros(1,nT), linspace(0,1,Nt-nT+1), ones(1,nT) ];
splineData.innerKnotsT = splineData.knotsT(nT+1:end-nT);

% Construct knots for end curve
splineData.endCurve.noPts = splineData.numVertices;
noPts = splineData.endCurve.noPts;

    if curveClosed % if end curve is closed

        splineData.endCurve.knotsEndCurve = linspace(0, 2*pi, noPts+1);
        splineData.endCurve.knotsEndCurve = splineData.endCurve.knotsEndCurve(1:end-1); % remove the last repeated point

    else % if end curve is open

        splineData.endCurve.knotsEndCurve = linspace(0, 2*pi, noPts);

    end

end

