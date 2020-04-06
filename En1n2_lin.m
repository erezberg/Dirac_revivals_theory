% function E = En1n2_lin(n,param)
% compute E(n), where n is a vector of densities of the four flavors, for
% linear DOS model
function E = En1n2_lin(n,param)
Ek = (2/3)*param.W*abs(n).^1.5;
E = param.U*(n'*param.utmat*n);
E = E + sum(Ek) - param.mu*sum(n);
if (param.B~=0)
    E = E - param.B*(param.spin'*n);
end
end
