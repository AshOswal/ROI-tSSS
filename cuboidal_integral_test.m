function cuboidal_integral_test

% This function creates a cuboid in spherical co-ordinates and then
% plots in cartesian co-ordinates. The cube is divided into 4 pyramids. 
% Each rectangle has a top, middle and bottom component as per the figure.
% The script also uses a quadrature based method to verify the volume (element)
% using integration in the spherical co-ordinate system. 

% CONVENTION NOTE. Phi is azimuth angle measured in radians. Theta is
% elevation from negative Z axis. This convention is carried forward in
% sph2cart

% Sx Sy Sz are the dimensions of the cuboid

Sx        = 10;
Sy        = 4;
Sz        = 2;

S         = [Sx Sy Sz];
Se        = 16;
fun       = @(phi,theta,r) r.^2.*sin(theta);

phi_min   = [-atan(Sy/Sx) atan(Sy/Sx) pi-atan(Sy/Sx) pi+atan(Sy/Sx)];
phi_max   = [atan(Sy/Sx) pi-atan(Sy/Sx) pi+atan(Sy/Sx) (2*pi)-atan(Sy/Sx)];
phi_mid   = [0     pi/2   pi     6/4*pi];
ratios    = [Sx/Sz Sy/Sz Sx/Sz Sy/Sz];

% note in this case that AA will vary
%% Middle

% ratio here should be Sx/Sz
% AA should be Sx
theta_min = @(phi,pm,ratio) atan(sec(phi-pm) *ratio);
theta_max = @(phi,pm,ratio) pi- atan(sec(phi-pm)*ratio);
r         = @ (phi,theta,AA,pm) AA./(2.*cos(phi-pm).*sin(theta));


%% Top

% AA should be Sz
theta_min_t = @(phi,pm,ratio) pi- atan(sec(phi-pm)*ratio);
theta_max_t = pi;
r_t         = @ (phi,theta,AA) AA./(2.*cos(pi-theta));

%% Bottom

% AA should be Sz
theta_min_b = 0;
theta_max_b = @(phi,pm,ratio) atan(sec(phi-pm) *ratio);
r_b         = @ (phi,theta,AA) AA./(2.*cos(theta));

%%

figure;

for n = 1:numel(phi_min)
phi         = repmat(linspace(phi_min(n),phi_max(n),Se),Se,1);
theta_MIN   = theta_min(phi(1,:),phi_mid(n),ratios(n));
theta_MAX   = theta_max(phi(1,:),phi_mid(n),ratios(n));
       
theta_MIN_t = theta_min_t(phi(1,:),phi_mid(n),ratios(n));
theta_MAX_b = theta_max_b(phi(1,:),phi_mid(n),ratios(n));

thetas    = []; thetas_b = []; thetas_t = [];
for j     = 1:numel(theta_MIN)
          thetas   = [thetas,linspace(theta_MIN(j),theta_MAX(j),Se)'];
          thetas_b = [thetas_b,linspace(theta_min_b,theta_MAX_b(j),Se)'];
          thetas_t = [thetas_t,linspace(theta_MIN_t(j),theta_max_t,Se)'];
end
%%
if ismember(n,[1,3])
    Mx      = v2mat(1:Sx,size(phi));
    xint    = Sx;
else
    Mx      = v2mat(1:Sy,size(phi));
    xint    = Sy;
end
Mz      = v2mat(1:Sz,size(phi));

phi_m   = repmat(phi,[1 1 size(Mx,3)]);
phi_t   = repmat(phi,[1 1 size(Mz,3)]);


thetas  = repmat(thetas,[1 1 size(Mx,3)]);
thetas_b  = repmat(thetas_b,[1 1 size(Mz,3)]);
thetas_t  = repmat(thetas_t,[1 1 size(Mz,3)]);

R         = r(phi_m,thetas,Mx,phi_mid(n));
R_t       = r_t(phi_t,thetas_t,Mz);
R_b       = r_b(phi_t,thetas_b,Mz);

[x,y,z]      = sph2cart(phi_m,thetas-pi/2,R);
[xt,yt,zt]   = sph2cart(phi_t,thetas_t-pi/2,R_t);
[xb,yb,zb]   = sph2cart(phi_t,thetas_b-pi/2,R_b);

for pvec = n:4
subplot(2,2,pvec);
plot3(x(:),y(:),z(:),'b.','MarkerFaceColor','b');hold on;
plot3(xt(:),yt(:),zt(:),'b.','MarkerFaceColor','b');hold on;
plot3(xb(:),yb(:),zb(:),'b.','MarkerFaceColor','b');hold on;
plot3(0,0,0,'g.','MarkerSize',15);
axis([-6 6 -4 4 -2 2]);
xlabel('X cm');ylabel('Y cm');zlabel('Z cm');
set(gca,'FontSize',12,'XGrid','on','YGrid','on','ZGrid','on')
end
% verifies that the volume is indeed 1/4th that of a cube

QTm(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n),ratios(n)),@(phi) theta_max(phi,phi_mid(n),ratios(n)),0,@ (phi,theta) r(phi,theta,xint,phi_mid(n)),'AbsTol',eps,'RelTol',1e-10);
QTt(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n),ratios(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,Sz),'AbsTol',eps,'RelTol',1e-10);
QTb(n) = integral3(fun,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n),ratios(n)),0,@ (phi,theta) r_b(phi,theta,Sz),'AbsTol',eps,'RelTol',1e-10);



end

S = sum(QTm+QTb+QTt);
disp(['Total calculated volume ' num2str(S)]);
disp(['Input volume ' num2str(prod(S))]);

function mat = v2mat(v,siz)
mat = [];
for n = 1:numel(v)
    mat = cat(numel(siz)+1,mat,repmat(v(n),siz));
end