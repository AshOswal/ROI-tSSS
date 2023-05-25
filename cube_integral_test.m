function cube_integral_test

% This function creates a cube based in spherical co-ordinates and then
% plots in cartesian co-ordinates. The cube is divided into 4 pyramid. 
% Each pyramid has a top, middle and bottom component as per the figure.
% The script also uses a quadrature based method to verify the volume (element)
% using integration in the spherical co-ordinate system. 

% CONVENTION NOTE. Phi is azimuth angle measured in radians. Theta is
% elevation from negative Z axis. This convention is carried forward in
% sph2cart

phi_min   = [-pi/4 pi/4   3/4*pi 5/4*pi];
phi_max   = [pi/4  3/4*pi 5/4*pi 7/4*pi];
phi_mid   = [0     pi/2   pi     6/4*pi];

A         = 5;
S         = 10;
fun       = @(phi,theta,r) r.^2.*sin(theta);

%% Middle
theta_min = @(phi,pm) atan(sec(phi-pm)); 
theta_max = @(phi,pm) pi-atan(sec(phi-pm));
r         = @ (phi,theta,AA,pm) AA./(2.*cos(phi-pm).*sin(theta));

%% Top
theta_min_t = @(phi,pm) pi - atan(sec(phi-pm));
theta_max_t = pi;
r_t         = @ (phi,theta,AA) AA./(2.*cos(pi-theta));

%% Bottom

theta_min_b = 0;
theta_max_b = @(phi,pm) atan(sec(phi-pm));
r_b         = @ (phi,theta,AA) AA./(2.*cos(theta));

%%

figure;

for n = 1:numel(phi_min)
phi         = repmat(linspace(phi_min(n),phi_max(n),S),S,1);
theta_MIN   = theta_min(phi(1,:),phi_mid(n));
theta_MAX   = theta_max(phi(1,:),phi_mid(n));
       
theta_MIN_t = theta_min_t(phi(1,:),phi_mid(n));
theta_MAX_b = theta_max_b(phi(1,:),phi_mid(n));

thetas    = []; thetas_b = []; thetas_t = [];
for j     = 1:numel(theta_MIN)
          thetas   = [thetas,linspace(theta_MIN(j),theta_MAX(j),S)'];
          thetas_b = [thetas_b,linspace(theta_min_b,theta_MAX_b(j),S)'];
          thetas_t = [thetas_t,linspace(theta_MIN_t(j),theta_max_t,S)'];
end
%%
M       = v2mat(1:5,size(phi));
phi     = repmat(phi,[1 1 size(M,3)]);
thetas  = repmat(thetas,[1 1 size(M,3)]);
thetas_b  = repmat(thetas_b,[1 1 size(M,3)]);
thetas_t  = repmat(thetas_t,[1 1 size(M,3)]);

R         = r(phi,thetas,M,phi_mid(n));
R_t       = r_t(phi,thetas_t,M);
R_b       = r_b(phi,thetas_b,M);

[x,y,z]      = sph2cart(phi,thetas-pi/2,R);
[xt,yt,zt]   = sph2cart(phi,thetas_t-pi/2,R_t);
[xb,yb,zb]   = sph2cart(phi,thetas_b-pi/2,R_b);

plot3(x(:),y(:),z(:),'ro','MarkerFaceColor','r');hold on;
plot3(xt(:),yt(:),zt(:),'bo','MarkerFaceColor','b');hold on;
plot3(xb(:),yb(:),zb(:),'ko','MarkerFaceColor','k');hold on;

% verifies that the volume is indeed 1/4th that of a cube

QTm(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A,phi_mid(n)),'AbsTol',eps,'RelTol',1e-10);    
QTt(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A),'AbsTol',eps,'RelTol',1e-10);    
QTb(n) = integral3(fun,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A),'AbsTol',eps,'RelTol',1e-10);    

end
S = sum(QTm+QTb+QTt);
disp(['Total volume ' num2str(S)]);

function mat = v2mat(v,siz)
mat = [];
for n = 1:numel(v)
    mat = cat(numel(siz)+1,mat,repmat(v(n),siz));
end