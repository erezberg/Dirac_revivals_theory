% function E = En1n2(n,param)
% compute E(n), where n is a vector of densities of the four flavors 
% (general DOS)
function E = En1n2(n,param)
Ek = interp1(param.n,param.Ek,n);
E = param.U*(n'*param.utmat*n);
E = E + sum(Ek) - param.mu*sum(n);
if (param.B~=0)
    E = E - param.B*(param.spin'*n);
end
end
