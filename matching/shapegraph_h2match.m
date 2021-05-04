% shapegraph_h2match.m: function for performing shape graph registration using second order elastic Sobolev metrics (H2 metrics).
%
% Input:
%   source          : structure containing source shape graph
%   target          : structure containing target shape graph
%   objfun          : structure containing parameters for the discretized matching objective function (optional)
%   optim           : structure containing parameters for the optimization procedure (optional)
%   initPath        : initialization for optimizing over path of shape graphs (optional, default is constant path of source)
%   initSourceRho   : initialization for optimizing over transformed source edge weights (optional, default is zero weight on all transformed source edges)
%   initTargetRho   : initialization for optimizing over transformed target edge weights (optional, default is zero weight on all transformed target edges)
%   initBeta        : initialization for optimizing over rotations (optional, default is zero rotation angle - handles curves in R^2 only)
%   initV           : initialization for optimizing over translations (optional, default is zero vector in R^d)
%   initKappa       : initialization for optimizing over scalings (optional, default is scaling factor 1)
%
% Output:
%   optPath         : structure containing optimal deformation path for transformed source (geodesic between source and target)
%   transfSource    : structure containing transformed source
%   transfTarget    : structure containing transformed target
%   updatedSource   : structure containing source with reordered vertices and component curves after applying Tarjan's algorithm
%   summary         : structure containing energies, costs and information from the optimization process


function [optPath, transfSource, transfTarget, updatedSource, summary] = shapegraph_h2match(source, target, varargin)


global template data gradavg


% Dimensions
[~, d] = size(source.x);


%------------------------------%
%      PRELIMINARY CHECKS      %
%------------------------------%

% Preliminary checks for well-defined vertices and connectivity matrix
if ~isfield(source, 'x')
    error('Provide source.x, i.e., vertex list for source shape graph.') 
end

if ~isfield(source, 'G')
    error('Provide source.G, i.e., edge list for source shape graph.') 
end

if ~isfield(target, 'x')
    error('Provide target.x, i.e., vertex list for target shape graph.') 
end

if ~isfield(target, 'G')
    error('Provide target.G, i.e., edge list for target shape graph.') 
end


%----------------------------------%
%      SET DEFAULT PARAMETERS      %
%----------------------------------%

% Set default signal, edge weights and topology for the source and target shape graphs
if ~isfield(source, 'f')
    source.f = zeros(size(source.x,1), 1);
end

if ~isfield(source, 'rho')
    source.rho = zeros(size(source.G,1), 1);
end

if ~isfield(source, 'topology')
    source.topology = 'general';
end

if ~isfield(target, 'f')
    target.f = zeros(size(target.x,1), 1);
end

if ~isfield(target, 'rho')
    target.rho = zeros(size(target.G,1), 1);
end

if ~isfield(target, 'topology')
    target.topology = 'general';
end

% Allocate structures for objective function and optimization procedure parameters
objfun = struct();
optim = struct();

% Use user-specified parameters
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii}, 'char'))
        switch (lower(varargin{ii}))
            case 'objfun'
                ii = ii + 1;
                objfun = varargin{ii};
            case 'optim'
                ii = ii + 1;
                optim = varargin{ii};
        end
    end
    ii = ii + 1;
end

% Add remaining default parameters in objfun and optim structures
[objfun, optim] = shapegraph_h2match_setdefaults(objfun, optim);


%--------------------------%
%     INITIALIZATIONS      %
%--------------------------%

% By default, initialize matching (with appropriate structure for calling Hanso bfgs code) by defining:
%   a) initial path of shape graphs (set to constant path [of each component curve] of source shape graph)
%   b) initial transformed source edge weights (set to source edge weights)
%   c) initial transformed target edge weights (set to target edge weights)
%   d) initial rotation angle (set to zero)
%   e) initial translation vector (set to zero vector)
%   f) initial scaling factor (set to 1)

init = struct();
init.path = [];
init.sourceRho = source.rho;
init.targetRho = target.rho;
init.beta = 0;
init.v = zeros(d,1);
init.kappa = 1;
         
% Update default initializations using user-specified initializations (if provided)
ii = 1;
while ii <= length(varargin)
    if (isa(varargin{ii}, 'char'))
        switch (lower(varargin{ii}))
            case 'initpath'
                ii = ii + 1;
                init.path = varargin{ii};
            case 'initsourcerho'
                ii = ii + 1;
                init.sourceRho = varargin{ii};
            case 'inittargetrho'
                ii = ii + 1;
                init.targetRho = varargin{ii};
            case 'initbeta'
                ii = ii + 1;
                init.beta = varargin{ii};
            case 'initv'
                ii = ii + 1;
                init.v = varargin{ii};
            case 'initkappa'
                ii = ii + 1;
                init.kappa = varargin{ii};  
            case {'objfun', 'optim'}
                ii = ii + 1;
            otherwise
                error('Invalid option: ''%s''.', varargin{ii});
        end  
    end
    ii = ii + 1;
end


%------------------------------------------------------------------------------------%
%     EXTRACT COMPONENTS OF SOURCE SHAPE GRAPH & CONSTRUCT PRODUCT TENSOR SPLINE     %
%------------------------------------------------------------------------------------%

% Split source shape graph into its components, i.e., into disjoint union of open curves (or single closed curve)
source_spline = source;
source_spline.rho = init.sourceRho;  % use user-specified source edge weights if provided
[connComp, updatedSource, cutVertices] = get_conn_comp(source_spline); 

% Build spline on each component of the source shape graph, store in cell array
splineDataNew = cell(length(connComp), 1);  % cell array of spline parameters for each component
curvePath = cell(length(connComp), 1);  % cell array of control points for deformation path of each component of transformed source shape graph
updatedSource.cPts = [];  % control points of source shape graph

for k = 1:length(connComp)
    
    % Define spline parameters for each component
    splineDataNew{k} = constructSplineData;
    splineDataNew{k}.dSpace = d;  % dimension of ambient space, R^d
    splineDataNew{k}.N = ceil(0.75*size(connComp{k}.x, 1));  % number of control points for spatial discretization
    splineDataNew{k}.numVertices = length(connComp{k}.x);  % number of vertices (to be used when interpolating transformed source)
    
    % Load default spline parameters if not pre-specified
    if ~isfield(objfun,'splineData')
        objfun.splineData = constructSplineData;
    end
    
    splineDataNew{k}.nS = objfun.splineData.nS;  % default nS = 2, quadratic splines are needed for H2 metrics
    splineDataNew{k}.Nt = objfun.splineData.Nt;  % default Nt = 10, number of control points for time discretization
    splineDataNew{k}.nT = objfun.splineData.nT;  % default Nt = 1, linear splines for time discretization 
    splineDataNew{k}.curveClosed = 1*strcmp(connComp{k}.topology, 'closed');
    
    % Add knot sequences and collocation matrices (quadData, quadDataTensor)
    splineDataNew{k} = finishSplineSetup(splineDataNew{k});
    
    % Construct spatial spline control points for each component curve
    connComp{k}.cPts = constructSplineApproximation(connComp{k}.x, splineDataNew{k});
    updatedSource.cPts = [updatedSource.cPts ; connComp{k}.cPts]; % update source structure
    
    % Construct product tensor spline (constant path) for each component curve of the source shape graph
    curvePath{k} = linearPath(connComp{k}.cPts, connComp{k}.cPts, splineDataNew{k});  
    
end

% Use default initialization for path of shape graphs if no user-specified path is provided
if isempty(init.path)
    init.path = curvePath;
end

% Update source structure with component curves of source shape graph (with user-specified source edge weights)
updatedSource.connComp = connComp;

% Construct gradient averaging matrix (with respect to control points, not vertices) for transformed source shape graph
[~, gradavg] = get_grad_avg(updatedSource, cutVertices);

% Update spline parameter structure (cell array of structs for each component of source shape graph)
objfun.splineData = splineDataNew;


%-------------------------%
%     ERROR MESSAGE       %
%-------------------------%

% Display error message if optimization over scalings is ill-defined
if optim.options.scaleInv == false && optim.options.optScal == true
    disp('optScal is only a valid option for scale invariant metrics.') 
       return
end


%--------------------------%
%    PERFORM MATCHING      %
%--------------------------%

% Set global variables
template = updatedSource; 
data = target;
data.rho = init.targetRho;  % use user-specified target edge weights if provided

% Perform matching
[optPath, transfSource, transfTarget, summary] = shapegraph_h2match_optimize(init, objfun, optim);

% Update source structure with component curves of source shape graph (with original source edge weights)
connComp_updatedSource = get_conn_comp(source);
updatedSource.connComp = connComp_updatedSource;


end

