%% setupQuadData
%
% Computes the collocation matrices for quadrature knots. Computes only
% those matrices for which information is set in splineData.
%
% Input
%   splineData
%       Information about splines used and quadrature degrees.
%
% Output
%   splineData
%       Adds quadData and quadDataTensor as fields in splineData.
%
% Note: This function is part of the H2 metrics library (https://github.com/h2metrics/h2metrics).


function splineData = setupQuadData( splineData )

% Set up quadData structure
quadData = struct('quadPointsS', [], 'quadPointsT', [], ...
    'noQuadPointsS', [], 'noQuadPointsT', [], ...
    'quadWeightsS', [], 'quadWeightsT', [], ...
    'B_S', [], 'Bu_S', [], 'Buu_S', [] ,'Buuu_S', [], ...
    'B_T', [], 'Bt_T', [],...
    'B_interpolS', [], 'B_endCurveS', []);

% Determine which collocation matrices need to be set up
curveClosed = splineData.curveClosed;
quadDegree = splineData.quadDegree;

doS = ~isempty(splineData.nS) && ~isempty(splineData.N) && ...
    ~isempty(splineData.quadDegree);
doT = ~isempty(splineData.nT) && ~isempty(splineData.Nt) && ...
    ~isempty(splineData.quadDegree) && length(splineData.quadDegree) >= 2;
doInterpolS = ~isempty(splineData.nS) && ~isempty(splineData.N) && ...
    ~isempty(splineData.interpolS);
doEndCurve = ~isempty(splineData.endCurve.noPts) && ~isempty(splineData.endCurve.knotsEndCurve);

% Compute collocation matrix for quadrature knots in space
if doS
    N = splineData.N;
    nS = splineData.nS;
    knotsS = splineData.knotsS;
    innerKnotsS = splineData.innerKnotsS;
    
    [quadPointsS, quadWeightsS] = gaussianQuadratureData( ...
        unique(innerKnotsS), 'degree', quadDegree(1) );
    noQuadPointsS = length(quadPointsS);
    
    quadData.quadPointsS = quadPointsS';
    quadData.noQuadPointsS = noQuadPointsS;
    quadData.quadWeightsS = quadWeightsS';
    
    noSder = min(4, nS+1);
    
    B_S_quad = spcol( knotsS, nS+1, ...
                      brk2knt( quadPointsS, noSder ), 'sparse');
    if curveClosed
        B_S_quad = [ B_S_quad(:,1:nS) + B_S_quad(:,end-nS+1:end), ...
                     B_S_quad(:,nS+1:end-nS) ];
    end
    
    quadData.B_S = B_S_quad(1:noSder:end, :);
    quadData.Bu_S = B_S_quad(2:noSder:end, :);
    if nS >= 2
        quadData.Buu_S = B_S_quad(3:noSder:end, :);
    else 
        quadData.Buu_S = zeros(size(B_S_quad(3:noSder:end,:)));
        quadData.Buu_S = sparse(quadData.Buu_S);
    end
    if nS >= 3
        quadData.Buuu_S = B_S_quad(4:noSder:end, :);
    end
end

% Compute collocation matrix for quadrature knots in time
if doT
    nT = splineData.nT;
    Nt = splineData.Nt;
    knotsT = splineData.knotsT;
    innerKnotsT = splineData.innerKnotsT;

    [quadPointsT, quadWeightsT] = gaussianQuadratureData( ...
        unique(innerKnotsT), 'degree', quadDegree(2) );

    noQuadPointsT = length(quadPointsT);
    quadData.quadPointsT = quadPointsT';
    quadData.noQuadPointsT = noQuadPointsT;
    quadData.quadWeightsT = quadWeightsT';

    noTder = 2;
    B_T_quad = spcol( knotsT, nT+1, ...
                      brk2knt( quadPointsT, noTder ), 'sparse');              
    quadData.B_T = B_T_quad(1:noTder:end,:);
    quadData.Bt_T = B_T_quad(2:noTder:end,:);
end

% Compute matrix for spline interpolation in space
if doInterpolS
    nS = splineData.nS;
    interpolS = splineData.interpolS;
    B_interpolS = spcol( splineData.knotsS, nS+1, ...
                         brk2knt(interpolS, 1), 'sparse');
    if curveClosed
        B_interpolS = [ B_interpolS(:,1:nS) + B_interpolS(:,end-nS+1:end), ...
                        B_interpolS(:,nS+1:end-nS) ];
    end
    quadData.B_interpolS = B_interpolS;
end

% Compute outer products for product tensor B-spline
quadDataTensor = struct('B',[],'Bu',[],'Buu',[],'Bt',[],'But',[],...
    'Buut',[],'Buuu',[], 'quadWeights',[],...
    'BuTr',[],'BuuTr',[],'BtTr',[],'ButTr',[],'BuutTr',[],...
    'B_phi',[],'Bu_phi',[],'Bt_phi',[],'Buu_phi',[],'But_phi',[],...
    'Buut_phi',[],'Buuu_phi',[],'nnz',[]) ;

if doS && doT
    quadDataTensor.quadWeights = reshape(quadWeightsS'*quadWeightsT,[],1);
    
    quadDataTensor.B = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 1, 1, splineData );
    quadDataTensor.Bu = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 2, 1, splineData );
    quadDataTensor.Bt = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 1, 2, splineData );
    quadDataTensor.But = createTensorCollocationMatrix( ...
        quadPointsS, quadPointsT, 2, 2, splineData );
    
    quadDataTensor.BuTr = quadDataTensor.Bu';
    quadDataTensor.BtTr = quadDataTensor.Bt';
    quadDataTensor.ButTr = quadDataTensor.But';
    quadDataTensor.nnz = nnz(quadDataTensor.Bu);
    
    if nS >= 2
        quadDataTensor.Buu = createTensorCollocationMatrix( ...
            quadPointsS, quadPointsT, 3, 1, splineData );
        
        quadDataTensor.Buut = createTensorCollocationMatrix( ...
            quadPointsS, quadPointsT, 3, 2, splineData );
        
        quadDataTensor.BuuTr = quadDataTensor.Buu';
        quadDataTensor.BuutTr = quadDataTensor.Buut';
    else 
        quadDataTensor.Buu = zeros(size(quadDataTensor.Bu));
        quadDataTensor.Buu = sparse(quadDataTensor.Buu);
        quadDataTensor.Buut = zeros(size(quadDataTensor.Bu));
        quadDataTensor.Buut = sparse(quadDataTensor.Buut);
        quadDataTensor.BuuTr = quadDataTensor.Buu';
        quadDataTensor.BuuTr = sparse(quadDataTensor.BuuTr);
        quadDataTensor.BuutTr = quadDataTensor.Buut';
        quadDataTensor.BuutTr = sparse(quadDataTensor.BuutTr);
    end
    if nS >= 3
        quadDataTensor.Buuu = createTensorCollocationMatrix( ...
            quadPointsS, quadPointsT, 4, 1, splineData );
    end
end

% Compute collocation matrix for end curve B-spline
if doEndCurve
    pts = splineData.endCurve.knotsEndCurve;
    B_endCurveS = spcol( splineData.knotsS, nS+1, ...
                    brk2knt(pts, 1), 'sparse'); 
                
    % Periodicity
    if curveClosed
        nS = splineData.nS;
        B_endCurveS = [B_endCurveS(:,1:nS) + B_endCurveS(:,end-nS+1:end),...
            B_endCurveS(:,nS+1:end-nS) ];
    end
           
    quadData.B_endCurveS = B_endCurveS;
    
end


% Save quadData and quadDataTensor as fields in splineData
splineData.quadData = quadData;
splineData.quadDataTensor = quadDataTensor;
    
end


