%% Set Directory
clear all; close all;
addpath(genpath('../'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ShapeGraph_H2match demo file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% source: curve with three branches
% target: straight curve with one straight branch
% model: h2 metric, weights on source only
% method: sfista


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Define Source & Target %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the source and target shape graphs:
%
% • The source shape graph should be defined as a structure with the following fields:
%       - source.x: a list of N vertices in R^d for the source shape graph [Nxd array] (mandatory)
%       - source.G: a list of M edges connecting the vertices [Mx2 array] (mandatory)
%       - source.rho: a list of edge weights [Mx1 array] (optional, default is zero edge weights)
%
% • Likewise, the target shape graph should be a structure with the fields target.x for
%   its vertices, target.G for its edges and optionally, target.rho for edge weights.

% Load data
load('../datasets/branches.mat')

% Source shape graph: curve with three branches
source = branch1;

% Target shape graph: straight curve with one straight branch
target = branch2;

% Plot source and target for viewing
figure()
hold on
plot_shape(target,'colorstyle','r.-','linewidth',1.5)  % target is red
plot_shape(source,'colorstyle','b.-','linewidth',1.5)  % source is blue
axis equal
hold off


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Select Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selecting parameters:
%
% • Select parameters for performing shape graph registration by defining the objfun and
%   optim structures. 
%  
% • The objfun structure contains parameters for the discretized relaxed
%   matching objective function, with the following fields:
%
%       - objfun.edge_weight: model selection parameter for varying edge weights. Possible
%         values are 'source', 'target', 'both' or 'fixed_weights'. (default is 'source')
%       - objfun.h2: structure containing coefficients for the H2 Riemannian metric.
%       - objfun.kernel_distance: structure containing varifold data attachment term parameters.
%       - objfun.tv: structure containing TV norm regularizer parameters.
%       - objfun.penalty: structure containing double-well penalty parameters.
%       - objfun.splineData: cell array of structures containing spline parameters for
%         discretizing the path of each component curve of the transformed source shape graph.
%
% • The optim structure contains parameters for the optimization procedure to minimize the
%   matching objective function, with the following fields:
%
%       - optim.method: optimization procedure parameter. Possible values are 'sfista' or 'bfgs'
%         (default is 'sfista')
%       - optim.options: structure containing parameters for optimization options.
%
% • Defining objfun and optim is optional, with default parameter values chosen by the 
%   shapegraph_h2match_setdefaults function when running the matching.
%
% • For more details on how to define objfun and optim, see the example right below.


% Model selection parameter
objfun.edge_weight = 'source';  % vary edge weights on 'source' (default), 'target', 'both' or 'fixed_weights'

% Optimization method
optim.method = 'sfista';  % 'sfista' (default) or 'bfgs'

% H2 energy parameters
objfun.h2.a0 = 0.1;  % L2 term (default: 0.1)
objfun.h2.a1 = 1;    % H1 term (default: 1)
objfun.h2.a2 = 1e-5; % H2 term (default: 0)
objfun.h2.a3 = 0;    % tangential H1 term (default: 0)
objfun.h2.a4 = 0;    % orthogonal H1 term (default: 0)

% Varifold parameters
%objfun.kernel_distance.distance = 'empty'; % varifold norm 
% ('empty' by default in which case we use 'varexpo'? 
% options are: 'var' (unoriented varifold with binet kernel) or 'cur' (current distance (oriented))
objfun.kernel_distance.kernel_geom = 'gaussian';  % functional form of geometric kernel (Possible values are? default: 'gaussian')
objfun.kernel_distance.kernel_signal = 'gaussian';   % what's this? (Possible values are? default is 'gaussian')
objfun.kernel_distance.kernel_grass = 'gaussian_oriented';  % functional form of grassmanian kernel (Possible values are 'gaussian_oriented', 'binet'? default: 'gaussian_oriented') 
objfun.kernel_distance.kernel_size_geom = .2;  % size of geometric kernel in data attachment term (default: 0.2*diam, where diam = 1)
objfun.kernel_distance.kernel_size_signal = 1;  % size of the functional kernel in the data attachment term (default: 1) 
objfun.kernel_distance.kernel_size_grass = 1.5;  % size of the Grassmanian kernel (default: 1.5)
objfun.kernel_distance.lambda_var = 3;  % weighting coefficient in front of data attachment term (default: 1)

% TV norm parameters
objfun.tv.lambda_tv_source = 0.02;  % weighting coefficient of TV norm for transformed source edge weights (default: 0.01)
objfun.tv.lambda_tv_target = 0.02;  % weighting coefficient of TV norm for target edge weights (default: 0.01)

% Double-well penalty parameters
objfun.penalty.lambda_pen_source = 0;  % weighting coefficient for double well penalty of transformed source weights (default: 1)
objfun.penalty.lambda_pen_target = 0;  % weighting coefficient for double well penalty of target weights (default: 1)
objfun.penalty.pen_min1 = -1;  % first minimizer of double well penalty (default: -1)
objfun.penalty.pen_min2 = 0;   % second minimizer of double well penalty (default: 0)

% Optimization options
optim.options.optGeom = true;  % optimize over geometric deformation of the transformed source (default: true)
optim.options.scaleInv = false;  % use scale invariant Riemannian metric (default: false)
optim.options.optRot = false;  % factor out rotations (default: false)
optim.options.optTra = false;  % factor out translations (default: false)
optim.options.optScal = false;  % factor out scalings (default: false)



%% %%%%%%%%%%%%%%%%%%%%%% Perform Shape Graph Matching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform shape graph matching by using the 'shapegraph_h2match' function.

% Run 1: Perform matching with the above parameter selection and default initializations
[optPath, transfSource, transfTarget, ~, summary1] = ...
                   shapegraph_h2match(source, target, 'objfun', objfun, 'optim', optim);

% Print costs and energies
summary1.costs

% Run 2: Perform matching with adjusted parameters and initializations from Run 1
objfun.kernel_distance.kernel_size_geom = 0.1; 
objfun.kernel_distance.lambda_var = 25; 
objfun.tv.lambda_tv_source = 0.1;
objfun.penalty.lambda_pen_source = 275;

[optPath, transfSource, transfTarget, updatedSource, summary] = ...
                 shapegraph_h2match(source, target, 'objfun', objfun, 'optim', optim,...
                                    'initPath', optPath,...
                                    'initSourceRho', transfSource.rho,...
                                    'initTargetRho', transfTarget.rho);

% Print costs and energies
summary.costs

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Visualize Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualize results from the matching process by using the 'shapegraph_h2match_viewresult' function.

close all
shapegraph_h2match_viewresult(optPath, transfSource, transfTarget, updatedSource, summary,...
          'save_fig', false, 'file_name', 'matching_branches_source_weights')

