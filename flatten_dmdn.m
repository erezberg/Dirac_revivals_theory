% [n_flat, dmdn_flat] = flatten_dmdn(n, mu, dndm)
% Take care of zeros of dmu/dn resulting from first order transitions. 
% This is done by detecting mu's where dn/dmu (computed by taking a numerical derivative)
% exceeds a certain threshold value. At these points, dn/dmu is set to 0. 
function [n_flat, dmdn_flat] = flatten_dmdn(n, mu, dndm)
    I = find(diff(n)./diff(mu)>5);
    n_flat = n;
    dmdn_flat = 1./dndm;
    for j = 1:length(I)
        n_flat = [n_flat(1:I(j)),n_flat(I(j)),n_flat(I(j)+1),n_flat((I(j)+1):end)];
        dmdn_flat = [dmdn_flat(1:I(j)),0,0,dmdn_flat(I(j)+1:end)];
        I = I + 2;
    end
end

