function Y = spharm_vec(theta,phi,l,m)
Y      = [];
thetas = theta(:);
phis   = phi(:);
for n  = 1:numel(thetas)
    y      = spharm(thetas(n),phis(n),l,m);
    Y      = [Y;y];
end
Y      = reshape(Y,size(theta));