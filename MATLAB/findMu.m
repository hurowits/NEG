% mu = findMu(s_vec,sigma)
% numerically solve for mu for every value in the vector s_vec    
function mu = findMu(s_vec,sigma)
for iS = 1:length(s_vec)
    s=s_vec(iS);
%     mu(iS) = fsolve(@(mu)s-1/mu*sinh(mu*sigma)/(mu*sigma),10);
    options = optimset('TolX',1e-11,'TolFun',1e-11);

   mu(iS) = fsolve(@(mu)s-1/mu*log(sinh(mu*sigma)/(mu*sigma)),10,options);
end
