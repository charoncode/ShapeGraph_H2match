% shapegraph_h2match_setdefaults.m: function to set default objective function and optimization procedure parameters for
% performing shape graph registration using second order elastic Sobolev metrics (H2 metrics).
%
% Input:
%   objfun : structure containing objective function parameters
%   optim  : structure containing optimization process parameters
%
% Output:
%   objfun and optim updated with default parameters

function [objfun, optim] = shapegraph_h2match_setdefaults(objfun, optim)


%% Computation method

comp_method = 'matlab';


%% Model selection parameters

if ~isfield(objfun, 'edge_weight') 
    objfun.edge_weight = 'source';  % vary edge weights on 'source', 'target', 'both' or 'fixed_weights'
end

if ~isfield(objfun, 'signal_type')
    objfun.signal_type = 'vertex';
end

if ~isfield(objfun, 'data_signal_type')
    objfun.data_signal_type = 'vertex';
end


%% H2 metric parameters

if ~isfield(objfun, 'h2')
    objfun.h2.a0 = 0.1; % L2 term
    objfun.h2.a1 = 1;   % H1 term      
    objfun.h2.a2 = 0;   % H2 term     
    objfun.h2.a3 = 0;   % tangential H1 term   
    objfun.h2.a4 = 0;   % orthogonal H1 term 
    
end

if isfield(objfun, 'h2') && ~isfield(objfun.h2, 'a0')
    objfun.h2.a0 = 0.1;
end

if isfield(objfun, 'h2') && ~isfield(objfun.h2, 'a1')
    objfun.h2.a1 = 1;
end

if isfield(objfun, 'h2') && ~isfield(objfun.h2, 'a2')
    objfun.h2.a2 = 0;
end

if isfield(objfun, 'h2') && ~isfield(objfun.h2, 'a3')
    objfun.h2.a3 = 0;
end

if isfield(objfun, 'h2') && ~isfield(objfun.h2, 'a4')
    objfun.h2.a4 = 0;
end


%% Varifold parameters

if ~isfield(objfun, 'distance')
    objfun.distance = 'kernel';
end

if ~isfield(objfun, 'kernel_distance')
    objfun.kernel_distance.kernel_geom = 'gaussian';
    objfun.kernel_distance.kernel_signal = 'gaussian';
    objfun.kernel_distance.kernel_grass = 'gaussian_oriented';  % options are 'gaussian_oriented', 'binet'
    objfun.kernel_distance.kernel_size_geom = .2;  % size of geometric kernel in data attachment term (default 0.2*diam, where diam = 1)
    objfun.kernel_distance.kernel_size_signal = 1;  % size of the functional kernel in the data attachment term (default 1) 
    objfun.kernel_distance.kernel_size_grass = 1.5;
    objfun.kernel_distance.method = comp_method;  % possible values are 'cuda' or 'matlab'
    objfun.kernel_distance.lambda_var = 1;  % weighting coefficient in front of data attachment term
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'kernel_geom')
    objfun.kernel_distance.kernel_geom = 'gaussian';  
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'kernel_signal')
    objfun.kernel_distance.kernel_signal = 'gaussian';  
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'kernel_grass')
    objfun.kernel_distance.kernel_grass = 'gaussian_oriented'; 
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'kernel_size_geom')
    objfun.kernel_distance.kernel_size_geom = .2; 
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'kernel_size_signal')
    objfun.kernel_distance.kernel_size_signal = 1; 
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'kernel_size_grass')
    objfun.kernel_distance.kernel_size_grass = 1.5; 
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'method')
    objfun.kernel_distance.method = comp_method; 
end
    
if isfield(objfun, 'kernel_distance') && ~isfield(objfun.kernel_distance, 'lambda_var')
    objfun.kernel_distance.lambda_var = 1; 
end


%% TV norm parameters

if ~isfield(objfun, 'tv')
    objfun.tv.lambda_tv_source = 0.01;  % weighting coefficient for tv norm for source weights
    objfun.tv.lambda_tv_target = 0.01;  % weighting coefficient for tv norm for target weights
end
    
if isfield(objfun, 'tv') && ~isfield(objfun.tv, 'lambda_tv_source')
    objfun.tv.lambda_tv_source = 0.01;
end
    
if isfield(objfun, 'tv') && ~isfield(objfun.tv, 'lambda_tv_target')
    objfun.tv.lambda_tv_target = 0.01;
end


%% Double-well penalty parameters

% Default is {0,1}-penalty
if ~isfield(objfun, 'penalty')
    objfun.penalty.lambda_pen_source = 1;  % weighting coefficient for double well penalty of transformed source weights
    objfun.penalty.lambda_pen_target = 1;  % weighting coefficient for double well penalty of target weights
    objfun.penalty.pen_min1 = -1;  % first minimizer of double well penalty 
    objfun.penalty.pen_min2 = 0;  % second minimizer of double well penalty
end
    
if isfield(objfun, 'penalty') && ~isfield(objfun.penalty, 'lambda_pen_source')
    objfun.penalty.lambda_pen_source = 1;
end
    
if isfield(objfun, 'penalty') && ~isfield(objfun.penalty, 'lambda_pen_target')
    objfun.penalty.lambda_pen_target = 1;
end
    
if isfield(objfun, 'penalty') && ~isfield(objfun.penalty, 'pen_min1')
    objfun.penalty.pen_min1 = -1;
end
    
if isfield(objfun, 'penalty') && ~isfield(objfun.penalty, 'pen_min2')
    objfun.penalty.pen_min2 = 0;    
end


%% Optimization parameters 

% Set default optimization method
if ~isfield(optim, 'method')
    optim.method = 'sfista';  % set sfista as default optimization method (alternatively set 'bfgs')
end

% Set default SFISTA parameters
if isfield(optim, 'method') && strcmp(optim.method, 'sfista')
    
    if ~isfield(optim, 'sfista')
        optim.sfista = struct();
    end
    
    if ~isfield(optim.sfista, 'prtlevel')
        optim.sfista.prtlevel = 1;  % 0 no printing, 1 prints last iteration, 2 prints at each iteration
    end
    
    if ~isfield(optim.sfista, 'maxit')
        optim.sfista.maxit = 2000;  % maximum number of iterations for each sfista iteration
    end
    
    if ~isfield(optim.sfista, 'maxit_gamma')
        optim.sfista.maxit_gamma = 25;  % maximum number of sfista iterations (i.e., number of times we increment gamma - the sfista gain parameter)
    end
    
    if ~isfield(optim.sfista, 'tol')
        optim.sfista.tol = 1e-6;  % stopping tolerance
    end
    
    if ~isfield(optim.sfista, 'gamma_source')
        optim.sfista.gamma_source = 1;  % sfista gain parameter for source weights
    end
    
    if ~isfield(optim.sfista, 'gamma_target')
        optim.sfista.gamma_target = 1;  % sfista gain parameter for target weights
    end
    
end
  
% Set default BFGS parameters (for minimization of the matching energy using smooth tv norm approximation)
if isfield(optim, 'method') && strcmp(optim.method, 'bfgs')
    
    if ~isfield(optim, 'bfgs')
        optim.bfgs = struct();
    end
    
    if ~isfield(optim.bfgs, 'prtlevel')
        optim.bfgs.prtlevel = 1;  % 0 no printing, 1 prints last iteration, 2 prints at each iteration
    end
    
    if ~isfield(optim.bfgs, 'maxit')
        optim.bfgs.maxit = 2000;  % maximum number of iterations for bfgs
    end
    
end
    
if ~strcmp(optim.method, 'sfista') && ~strcmp(optim.method, 'bfgs')
    disp('Invalid option for optim.method -- use ''sfista'' or ''bfgs''.')
end

    
%% Optimization options

if ~isfield(optim, 'options')
    optim.options.optGeom = 1;  % turn on geometric optimization by default
    optim.options.scaleInv = 0;  % use constant coefficient Riemannian metric by default
    optim.options.optRot = 0;  % no optimization over rotations by default
    optim.options.optTra = 0;  % no optimization over translations by default
    optim.options.optScal = 0;  % no optimization over scalings by default
end

if isfield(optim, 'options') && ~isfield(optim.options,'geom_opt')
    optim.options.optGeom = 1;
end
    
if isfield(optim, 'options') && ~isfield(optim.options, 'scaleInv')
    optim.options.scaleInv = 0;
end
    
if isfield(optim, 'options') && ~isfield(optim.options, 'optRot')
    optim.options.optRot = 0;
end
    
if isfield(optim, 'options') && ~isfield(optim.options, 'optTra')
    optim.options.optTra = 0;
end
    
if isfield(optim, 'options') && ~isfield(optim.options, 'optScal')
    optim.options.optScal = 0;
end

if ~isfield(optim, 'kernel_size_geom_reg')
    optim.kernel_size_geom_reg = 0;
end


end