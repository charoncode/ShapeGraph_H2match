% h2_energyPath.m: function to compute the H2 energy of a path of curves and its differential.
%
% Input: 
%   curvePath     : path of curves
%   objfun        : structure containing parameters for the discretized matching functional 
%   optim         : structure containing parameters for the optimization procedure
%   splineData    : structure containing spline parameters
%
% Output: 
%  energyPath     : H2 energy of path of curves
%  dENR_curvePath : derivative of H2 energy wrt to the path of curves, except the source at t=0 and the endcurve at t=1,
%                   i.e., wrt to [c_2,...,c_{Nt-1}]
%  dENR_endCurve  : derivative of the energy of path of curves wrt to the endcurve c(1)/c_{Nt}


function [energyPath, dENR_curvePath, dENR_endCurve] = h2_energyPath(curvePath, objfun, optim, splineData)

% Dimensions
d = size(curvePath, 2);

% H2 metric parameters
a0 = objfun.h2.a0; 
a1 = objfun.h2.a1; 
a2 = objfun.h2.a2; 
a1v = objfun.h2.a3; 
a1n = objfun.h2.a4;

% Optimization parameters
options = optim.options;
scaleInv = options.scaleInv;

% Spline parameters
N = splineData.N; 
Nt = splineData.Nt;
quadDataTensor = splineData.quadDataTensor;
noQuadPointsS = splineData.quadData.noQuadPointsS;
noQuadPointsT = splineData.quadData.noQuadPointsT;
quadWeightsS = splineData.quadData.quadWeightsS;

% Evaluate path and its derivatives at quadrature sites
Cu = quadDataTensor.Bu*curvePath;
Ct = quadDataTensor.Bt*curvePath;
Cut = quadDataTensor.But*curvePath;
Cuu = quadDataTensor.Buu*curvePath;
Cuut = quadDataTensor.Buut*curvePath;

CuCuu = sum(Cu.*Cuu,2);
CuCuu2 = CuCuu.^2;
CtCt = sum(Ct.*Ct,2);
CutCut = sum(Cut.*Cut,2);
CutCuut = sum(Cut.*Cuut,2);
CuutCuut = sum(Cuut.*Cuut,2);
CutCu = sum(Cut.*Cu,2);
CutCu2 = CutCu.^2;

Cspeed = sum( Cu.^2 , 2).^(1/2);

CspeedInv = 1./Cspeed;
CspeedInv2 = CspeedInv .* CspeedInv;
CspeedInv3 = CspeedInv2 .* CspeedInv;
CspeedInv5 = CspeedInv3 .* CspeedInv2;
CspeedInv7 = CspeedInv5 .* CspeedInv2;
CspeedInv9 = CspeedInv7 .* CspeedInv2;

% L2 energy terms
Ct_L2 = CtCt.*Cspeed;

% H1 energy terms
Ct_H1 = CspeedInv.*CutCut;
Ct_H1v = CutCu2.*CspeedInv3;

% H2 energy terms
Ct_H2 = CutCut .* CuCuu2 .* CspeedInv7 ...
    - 2 * CutCuut .* CuCuu .* CspeedInv5 ...
    + CuutCuut .* CspeedInv3;

if scaleInv == 0 % constant coefficient metric
    Ct_L2H1H2 = a0*Ct_L2 +(a1+a1n)*Ct_H1 + (a1v-a1n)*Ct_H1v + a2*Ct_H2;
    energyIntegrand = Ct_L2H1H2;
    
else
    % Compute length-weighted (scal-invariant) metric
    ellVec = quadWeightsS'*reshape(Cspeed,noQuadPointsS,noQuadPointsT);
    ellInvVec = ellVec.^(-1);
    ellInv3Vec = ellVec.^(-3);
    
    ell = reshape( repmat( ellVec , noQuadPointsS, 1), [],1) ;
    ellInv = reshape( repmat( ellInvVec , noQuadPointsS, 1), [],1);
    ellInv3 = reshape( repmat( ellInv3Vec , noQuadPointsS, 1), [],1);
    
    Ct_L2H1H2 = a0*ellInv3.*Ct_L2 +(a1+a1n)*ellInv.*Ct_H1 ...
        + (a1v-a1n)*ellInv.*Ct_H1v + a2*ell.*Ct_H2;
    energyIntegrand = Ct_L2H1H2;
    
end

% Compute final energy of path of curves
energyPath = quadDataTensor.quadWeights' * energyIntegrand;

% Compute differential of energy of path of curves wrt to spline control points of path of curves
if scaleInv == 0 % constant H2 coefficients
        term1 = a0*CspeedInv.*CtCt ...
            - (a1+a1n)*CspeedInv3.*CutCut ...
            - a2*7*CspeedInv9.*(CuCuu2).*CutCut ...
            + a2*10*CspeedInv7.*CuCuu.*CutCuut ...
            - a2*3*CspeedInv5.*CuutCuut ...
            - (a1v-a1n)*3*CspeedInv5.*CutCu2;
        term2 = a2*2*CspeedInv7.*CuCuu.*CutCut ...
            -a2*2*CspeedInv5.*CutCuut;
        term3 = 2*a0*Cspeed;
        term4 = 2*(a1+a1n)*CspeedInv + 2*a2*CspeedInv7.*CuCuu.^2;
        term5 = -a2*2*CspeedInv5.*CuCuu;
        term6 = 2*a2*CspeedInv3;
        term7 = 2*(a1v-a1n)*CspeedInv3.*CutCu;
        
        term1 = term1.*quadDataTensor.quadWeights;
        term2 = term2.*quadDataTensor.quadWeights;
        term3 = term3.*quadDataTensor.quadWeights;
        term4 = term4.*quadDataTensor.quadWeights;
        term5 = term5.*quadDataTensor.quadWeights;
        term6 = term6.*quadDataTensor.quadWeights;
        term7 = term7.*quadDataTensor.quadWeights;
        
        dENR_curvePath = zeros(N*Nt,d);
        
        for kk = splineData.dSpace:-1:1
            termBu = (term1.*Cu(:,kk) + term2.*Cuu(:,kk) + term7.*Cut(:,kk))'*quadDataTensor.Bu;
            termBuu = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
            termBt = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
            termBut = (term4.*Cut(:,kk)+term5.*Cuut(:,kk)+term7.*Cu(:,kk))'*quadDataTensor.But; %last term here is new
            termBuut = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;
            
            dENR_curvePath(:,kk) = termBu + termBuu + termBt + termBut + termBuut;
        end
  
    elseif scaleInv == 1 % scale invariant H2 metric
        ellInv2Vec = ellVec.^(-2);
        ellInv2 = reshape( repmat( ellInv2Vec , noQuadPointsS, 1), [],1);
        ellInv4Vec = ellVec.^(-4);
        ellInv4 = reshape( repmat( ellInv4Vec , noQuadPointsS, 1), [],1);
        
        Ct_L2_IntTheta = quadWeightsS'*reshape(Ct_L2,noQuadPointsS,noQuadPointsT);
        Ct_H1_IntTheta = quadWeightsS'*reshape(Ct_H1,noQuadPointsS,noQuadPointsT);
        Ct_H1v_IntTheta = quadWeightsS'*reshape(Ct_H1v,noQuadPointsS,noQuadPointsT);
        Ct_H2_IntTheta = quadWeightsS'*reshape(Ct_H2,noQuadPointsS,noQuadPointsT);
        
        Ct_L2_IntTheta = reshape( repmat( Ct_L2_IntTheta , noQuadPointsS, 1), [],1) ;
        Ct_H1_IntTheta = reshape( repmat( Ct_H1_IntTheta , noQuadPointsS, 1), [],1) ;
        Ct_H1v_IntTheta = reshape( repmat( Ct_H1v_IntTheta , noQuadPointsS, 1), [],1) ;
        Ct_H2_IntTheta = reshape( repmat( Ct_H2_IntTheta , noQuadPointsS, 1), [],1) ;
        
        term1 = a0*ellInv3.*CspeedInv.*CtCt ...
            - (a1+a1n)*ellInv.*CspeedInv3.*CutCut ...
            - a2*7*ell.*CspeedInv9.*(CuCuu2).*CutCut ...
            + a2*10*ell.*CspeedInv7.*CuCuu.*CutCuut ...
            - a2*3*ell.*CspeedInv5.*CuutCuut ...
            - (a1v-a1n)*3*ellInv.*CspeedInv5.*CutCu2 ...
            - a0*3*ellInv4.*Ct_L2_IntTheta.*CspeedInv ... %L-W L2 term
            - (a1+a1n)*ellInv2.*Ct_H1_IntTheta.*CspeedInv ...
            - (a1v-a1n)*ellInv2.*Ct_H1v_IntTheta.*CspeedInv ...
            + a2*Ct_H2_IntTheta.*CspeedInv;
        term2 = a2*2*ell.*CspeedInv7.*CuCuu.*CutCut ...
            -a2*2*ell.*CspeedInv5.*CutCuut;
        term3 = 2*a0*ellInv3.*Cspeed;
        term4 = 2*(a1+a1n)*ellInv.*CspeedInv + 2*a2*ell.*CspeedInv7.*CuCuu.^2;
        term5 = -a2*2*ell.*CspeedInv5.*CuCuu;
        term6 = 2*a2*ell.*CspeedInv3;
        term7 = 2*(a1v-a1n)*ellInv.*CspeedInv3.*CutCu;
        
        term1 = term1.*quadDataTensor.quadWeights;
        term2 = term2.*quadDataTensor.quadWeights;
        term3 = term3.*quadDataTensor.quadWeights;
        term4 = term4.*quadDataTensor.quadWeights;
        term5 = term5.*quadDataTensor.quadWeights;
        term6 = term6.*quadDataTensor.quadWeights;
        term7 = term7.*quadDataTensor.quadWeights;
        
        dENR_curvePath = zeros(N*Nt,d);
        
        for kk = d:-1:1
            termBu = (term1.*Cu(:,kk) + term2.*Cuu(:,kk) + term7.*Cut(:,kk))'*quadDataTensor.Bu;
            termBuu = (term2.*Cu(:,kk))'*quadDataTensor.Buu;
            termBt = (term3.*Ct(:,kk))'*quadDataTensor.Bt;
            termBut = (term4.*Cut(:,kk)+term5.*Cuut(:,kk)+term7.*Cu(:,kk))'*quadDataTensor.But;
            termBuut = (term6.*Cuut(:,kk)+term5.*Cut(:,kk))'*quadDataTensor.Buut;
            
            dENR_curvePath(:,kk) = termBu + termBuu + termBt + termBut + termBuut;
        end
end

% Differential of H2 energy wrt endcurve
dENR_endCurve = dENR_curvePath(end-N+1:end,:);

% Differential of H2 energy wrt to all other curves (except the source c0 at t=0 which is fixed)
dENR_curvePath = dENR_curvePath(N+1:end-N,:); 

end


