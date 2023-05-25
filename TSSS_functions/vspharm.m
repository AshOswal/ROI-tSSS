% X = vspharm(theta,phi,l,m)
%
% Vector spherical harmonic funtion (Arfken, p. 708) for quantum
% numbers l and m. Angles theta and phi should be given in radians.
% First component of X is the e_theta component and the second
% component is the e_phi component.
%
function X = vspharm(theta,phi,l,m)

scale = -1/sqrt(l*(l+1));
scale_sph = (-1)^m*sqrt((2*l+1)*prod(1:(l-m))/(4*pi*prod(1:(l+m))));
scale_minus = 1/((-1)^(m-1)*sqrt((2*l+1)*prod(1:(l-(m-1)))/(4*pi*prod(1:(l+(m-1))))));
scale_plus = 1/((-1)^(m+1)*sqrt((2*l+1)*prod(1:(l-(m+1)))/(4*pi*prod(1:(l+(m+1))))));

Y = spharm(theta,phi,l,m);
if m > -l
   Yminus = spharm(theta,phi,l,m-1);
else
   Yminus = 0;
end
if m < l
   Yplus = spharm(theta,phi,l,m+1);
else
   Yplus = 0;
end

dY = 0.5*scale_sph*((l+m)*(l-m+1)*scale_minus*Yminus*complex(cos(phi),sin(phi))-scale_plus*Yplus*complex(cos(-phi),sin(-phi)));  % dY/dtheta
if theta == 0 % To prevent division by zero
   theta = eps;
end
X(1) = (scale*m/sin(theta))*Y;
X(2) = scale*i*dY;

