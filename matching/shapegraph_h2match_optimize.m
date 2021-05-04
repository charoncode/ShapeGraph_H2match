% shapegraph_h2match_optimize.m: function for specifying optimization procedure to be used when performing shape graph 
% registration using second order elastic Sobolev metrics (H2 metrics). 
%
% Input:
%   init          : structure containing initializations for path of shape graphs, transformed source and/or
%                   transformed target edge weights, rotations, translations and scalings
%   objfun        : structure containing parameters for the discretized matching objective function
%   optim         : structure containing parameters for the optimization procedure
%
% Output:
%   optPath       : structure containing optimal deformation path of transformed source
%   transfSource  : structure containing transformed source
%   transfTarget  : structure containing transformed target
%   summary       : structure containing energies, costs and information from the optimization process


function [optPath, transfSource, transfTarget, summary] = shapegraph_h2match_optimize(init, objfun, optim)


global template data objfunc optimc


%-----------------------%
%      SET OPTIONS      %
%-----------------------%

list_of_variables = {'curvePath', 'rho0', 'rho1', 'beta', 'v', 'kappa'};

objfun = set_objfun_option(objfun, init, data);

switch lower(optim.method)
    
    case 'bfgs'  % use smooth TV norm approximation, optimize via BFGS
      optimc = set_optim_option(optim, objfun, list_of_variables, 'enr_shapegraph_h2match_bfgs');
      template.E_adj = edge_adjacency(template.G);  % edge adjacency matrix of source
      data.E_adj = edge_adjacency(data.G);  % edge adjacency matrix of target
      
    case 'sfista'  % optimize via SFISTA
      optimc = set_optim_option(optim, objfun, list_of_variables, 'enr_shapegraph_h2match_sfista'); 
      template.D = diff_operator(template);  % difference operator for computing TV norm of source edge weights
      data.D = diff_operator(data);  % difference operator for computing TV norm of target edge weights
      
end


%--------------------%
%      OPTIMIZE      %
%--------------------%

objfunc = objfun;
splineData = objfun.splineData;
initPath = [];

% Concatenate control points from the path of curves corresponding to each component of the source shape graph
for k = 1:length(template.connComp)
            
    Nk = splineData{k}.N;
    compPath = init.path{k};
    initPath = [initPath ; compPath(Nk+1:end, :)];  % all control points (except for the source at time t=0)
        
end

% Organize initializations into appropriate structure for calling Hanso BFGS
switch lower(objfunc.edge_weight)
    
    case 'both' % if optimizing over both source and target edge weights
        
        % organize variables into long (d*N*(Nt-1)+E0+E1+d+2)x1 array
        initvec = [initPath(:); ...
                   init.sourceRho(:); init.targetRho(:); ...
                   init.beta; init.v(:); init.kappa];
    
    case 'source' % if optimizing over source edge weights only
        
        % organize variables into long (d*N*(Nt-1)+E0+d+2)x1 array
        initvec = [initPath(:); ...
                   init.sourceRho(:); ...
                   init.beta; init.v(:); init.kappa];
    
    case 'target' % if optimizing over target edge weights only
        
        % organize variables into long (d*N*(Nt-1)+E1+d+2)x1 array
        initvec = [initPath(:); ...
                   init.targetRho(:);
                   init.beta; init.v(:); init.kappa];      

          
    case 'fixed_weights' % if not optimizing over weights (fixed weights)
        
        % organize variables into long (d*N*(Nt-1)+d+2)x1 array
        initvec = [initPath(:); ...
                   init.beta; init.v(:); init.kappa];           
        
end

% Choose optimization method and minimize matching objective function
switch lower(optim.method)
          
    case 'bfgs' % use BFGS with smooth TV norm approximation, return costs and transformed shape graphs
        [X, summary] = perform_bfgs(optimc.bfgs.enr_name, initvec, optimc.bfgs);
        [summary.costs, optPath, transfSource, transfTarget] = costs_shapegraph_h2match_bfgs(X);
           
    case 'sfista' % use SFISTA, return costs and transformed shape graphs
        [X,summary] = perform_h2sfista(initvec, objfunc, optimc);
        [summary.costs, optPath, transfSource, transfTarget] = costs_shapegraph_h2match_sfista(X);
                
end


%-------------------------%
%     SAVE PARAMETERS     %
%-------------------------% 

summary.parameters.objfun = objfun;
summary.parameters.optim = optim;


end

