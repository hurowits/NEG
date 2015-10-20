function pl = polylogN(n,x)

for iX=1:length(x)
    pl(iX)=double(feval(symengine, 'polylog', 2, x(iX)));
end