%% Set Directory
clear all; close all;
addpath(genpath('../'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ShapeGraph_H2match demo file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% source: swedish leaf with stem
% target: swedish leaf without stem and tip
% model: h2 metric, weights on source only
% method: sfista


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Define Source & Target %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
load('../datasets/SwedishLeaves.mat')

% source shape: swedish leaf with stem
N0 = 100;
source.x = ReSampleCurve(leafstem3',N0+1)'; % resample curve
source.x = source.x(2:end,:); % remove last repeated point
source.x = source.x - mean(source.x); % center at origin
diam_t = max(max(source.x) - min(source.x)); % calculate diameter
source.x = (source.x / diam_t) * 1.5; % scale diameter
source.G = [1:length(source.x) ; 2:length(source.x) 1]'; % define edge list (this is a closed curve)

% target shape: swedish leaf without stem and tip
N1 = 100;
target.x = ReSampleCurve(leaf3',N1+1)'; % resample curve
target.x = target.x(2:end,:); % remove last repeated point
target.x = target.x - mean(target.x); % center at origin
diam_s = max(max(target.x) - min(target.x)); % calculate diameter
target.x = target.x / diam_s ; % scale to unit diameter 
ind1 = 70:90; ind2 = 33:47;
target = remove_segment(target,{ind1, ind2}); % remove segment and define edge list

% Plot source and target for viewing
figure()
hold on
plot_shape(target,'colorstyle','r.-','linewidth',1.5)  % target is red
plot_shape(source,'colorstyle','b--','linewidth',1.5)  % source is blue
axis equal
hold off


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Select Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model selection parameter
objfun.edge_weight = 'source';  % vary edge weights on 'source', 'target', 'both' or 'fixed_weights'

% Optimization method
optim.method = 'sfista';  % 'bfgs' or 'sfista'

% H2 energy parameters
objfun.h2.a0 = 0.1;    % L2 term
objfun.h2.a1 = 1;      % H1 term      
objfun.h2.a2 = 1e-5;   % H2 term     
objfun.h2.a3 = 0;      % tangential H1 term   
objfun.h2.a4 = 0;      % orthogonal H1 term 

% Varifold parameters
objfun.kernel_distance.distance = 'var'; % varifold norm
objfun.kernel_distance.kernel_size_geom = 0.15; 
objfun.kernel_distance.kernel_size_grass = 1.5;
objfun.kernel_distance.lambda_var = 1;

% TV norm and penalty parameters
objfun.tv.lambda_tv_source = 0.01;
objfun.penalty.lambda_pen_source = 0;


%% %%%%%%%%%%%%%%%%%%%%%% Perform Shape Graph Matching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run 1: Perform matching with default parameter selections and default initializations
[optPath, transfSource, transfTarget, ~, summary1] = ...
                  shapegraph_h2match(source, target, 'objfun', objfun, 'optim', optim);

% Print costs and energies
summary1.costs

% Run 2: Perform matching with adjusted parameters and initializations from Run 1
objfun.kernel_distance.kernel_size_geom = .025; 
objfun.kernel_distance.lambda_var = 10;
objfun.tv.lambda_tv_source = 0;
objfun.penalty.lambda_pen_source = 250;

[optPath, transfSource, transfTarget, updatedSource, summary] = ...
                shapegraph_h2match(source, target, 'objfun', objfun, 'optim', optim,...
                                   'initPath', optPath,...
                                   'initSourceRho', transfSource.rho,...
                                   'initTargetRho', transfTarget.rho);

% Print costs and energies
summary.costs


%% %%%%%%%%%%%%%%%%%%%%%%%%%% Visualize Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
shapegraph_h2match_viewresult(optPath, transfSource, transfTarget, updatedSource,...
          summary, 'save_fig', false, 'file_name', 'matching_leaves_source_weights')