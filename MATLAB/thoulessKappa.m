% [kappa,signKappa,expKappa,E_vec] = thoulessKappa(E_j,n,g,E_vec)
% Calculates inverse localization length (2D electrostatic potential).
% Input: 
%   E_j - vector of charge locations
%   n - index of last charge to include in calculation
%   g - vector of coupling coefficients as defined in NEG 
%   E_vec - vector of values at which kappa is calculated
% Output:
%   All output is evaluated at points given by E_vec
%   kappa - inverse localization length evaluated at E_vec
%   signKappa - sign of the spectral determinant 
%   expKappa - spectral determinant evaluated at E_vec
%   E_vec is E_vec (used for debugging)
function [kappa,signKappa,expKappa,E_vec] = thoulessKappa(E_j,n,g,E_vec)

n=min(n,length(E_j));
% E_vec = [linspace(E_j(1),E_j(n)*1.1,1e5)]';
% E_vec = [linspace(E_j(470),E_j(500),1e5)]';

for iE=1:length(E_vec)
    E=E_vec(iE);
    kappa(iE) = sum(log(abs(E-E_j)))-sum(log(g));
    
    expKappa(iE)=prod(E-E_j);
    
    signKappa(iE) = prod(sign(E-E_j));
    
end
expKappa=expKappa/prod(g);