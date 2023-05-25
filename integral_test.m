clear all;
L  = [8 8];
M  = [5 5];
A  = 2;
%%%%%%%%%%%%%%% TEST FOR ORTHONORMALITY OF Spherical Harmonic Functions %%%%%%%%%%%%%%

% spharm_vec.m is vectorised wrapper for spharm.m

%fun = @(theta,phi) spharm_vec(theta,phi,L(1),M(1)).*conj(spharm_vec(theta,phi,L(2),M(2))).*sin(theta);
%Q = integral2(fun,0,pi,0,2*pi,'Method','iterated','AbsTol',1e-10,'RelTol',1e-10);
% Q yields unity as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fun   = @(phi,theta,r) sum(vsh_modified_in_AO(theta,phi,L(1),M(1)).*conj(vsh_modified_in_AO(theta,phi,L(2),M(2))),3) .*sin(theta).*r.^(sum(L)+2);
fun_t = @(phi,theta,r) r.^2.*sin(theta);
%% split into pyramid and upper and lower portions to make rectangle

phi_min   = [pi*7/4 pi/4   3/4*pi 5/4*pi];
phi_max   = [pi/4  3/4*pi 5/4*pi 7/4*pi];
phi_mid   = [0     pi/2   pi     6/4*pi];

% Middle
theta_min = @(phi,pm) atan(sec(phi-pm)); 
theta_max = @(phi,pm) pi-atan(sec(phi-pm));
r         = @ (phi,theta,AA,pm) AA./(2.*cos(phi-pm).*sin(theta));

% Top
theta_min_t = @(phi,pm) pi - atan(sec(phi-pm));
theta_max_t = pi;
r_t         = @ (phi,theta,AA) AA./(2.*cos(pi-theta));

% Bottom
theta_min_b = 0;
theta_max_b = @(phi,pm) atan(sec(phi-pm));
r_b         = @ (phi,theta,AA) AA./(2.*cos(theta));

for n = 1:numel(phi_min)
    
QTm(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A,phi_mid(n)));    
QTt(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A));    
QTb(n) = integral3(fun,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A));    

Q1(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A,phi_mid(n)));    
Q2(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A));    
Q3(n) = integral3(fun_t,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A));    

end

Y = sum(QTm+QTt+QTb);




