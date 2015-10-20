% W = generateRateMatrix(x,g)
% Input: 
%   x - stochastic field vector
%   g - couplings vector
% Output:
%   W - rate matrix 

function W = generateRateMatrix(x,g)

N = length(g);
upperDiag = g(1:N-1).*exp(x(1:N-1)/2);
lowerDiag = g(1:N-1).*exp(-x(1:N-1)/2);

W = diag(upperDiag,1)+diag(lowerDiag,-1);
W(1,N) = g(N)*exp(-x(N)/2);
W(N,1) = g(N)*exp(x(N)/2);

diagonal = - sum(W,1);

W = W + diag(diagonal);

