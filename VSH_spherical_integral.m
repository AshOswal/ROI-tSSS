function [Y,Y1] = VSH_spherical_integral(L,M,A)

% L is the L value of the matrix
% M is the M value of the matrix
% A is the side of the cubic ROI - in metres

if isempty(L)   
    L  = [8 8];
end
if isempty(M)
    M  = [5 5];
end
if isempty(A)
    A  = 0.1;
end
useV   = 0;
useX   = 1;
%%%%%%%%%%%%%%% TEST FOR ORTHONORMALITY OF (Vector) Spherical Harmonic Functions %%%%%%%%%%%%%%

% spharm_vec.m is vectorised wrapper for spharm.m
%{
fun = @(theta,phi) spharm_vec(theta,phi,L(1),M(1)).*conj(spharm_vec(theta,phi,L(2),M(2))).*sin(theta);
Q = integral2(fun,0,pi,0,2*pi,'Method','iterated','AbsTol',1e-10,'RelTol',1e-10);

fun = @(theta,phi) sum(vsh_modified_in_AO(theta,phi,L(1),M(1)).*conj(vsh_modified_in_AO(theta,phi,L(2),M(2))),3).*sin(theta);
Q1 = integral2(fun,0,pi,0,2*pi)

fun = @(theta,phi) sum(vspharm_AO(theta,phi,L(1),M(1)).*conj(vspharm_AO(theta,phi,L(2),M(2))),3).*sin(theta);
Q2 = integral2(fun,0,pi,0,2*pi)
%}
% Q yields unity as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if useV % uses Vlm in taulu paper
    fun   = @(phi,theta,r) sum(vsh_modified_in_AO(theta,phi,L(1),M(1)).*conj(vsh_modified_in_AO(theta,phi,L(2),M(2))),3) .*sin(theta).*(r.^(sum(L)+2));
elseif useX
    fun   = @(phi,theta,r) sum(conj(vspharm_AO(theta,phi,L(1),M(1))).* vspharm_AO(theta,phi,L(2),M(2)),3) .*sin(theta).*(r.^(sum(L)+2)); 
end
%fun_t = @(phi,theta,r) r.^2.*sin(theta);
%% split into pyramid and upper and lower portions to make rectangle

% phi_min   = [-pi/4 pi/4   3/4*pi 5/4*pi];
% phi_max   = [pi/4  3/4*pi 5/4*pi 7/4*pi];
% phi_mid   = [0     pi/2   pi     6/4*pi];

% Middle
phi_min   = 0;
phi_max   = 2*pi;
theta_min = 0;%@(phi,pm) atan(sec(phi-pm)); 
theta_max = pi;%@(phi,pm) pi-atan(sec(phi-pm));
r         = A;%@ (phi,theta,AA,pm) AA./(2.*cos(phi-pm).*sin(theta));

% % Top
% theta_min_t = @(phi,pm) pi - atan(sec(phi-pm));
% theta_max_t = pi;
% r_t         = @ (phi,theta,AA) AA./(2.*cos(pi-theta));
% 
% % Bottom
% theta_min_b = 0;
% theta_max_b = @(phi,pm) atan(sec(phi-pm));
% r_b         = @ (phi,theta,AA) AA./(2.*cos(theta));

%for n = 1:numel(phi_min)
    
% QTm(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A,phi_mid(n)),'AbsTol', 1e-6,'RelTol',1e-6);    
% QTt(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A),'AbsTol', 1e-6,'RelTol',1e-6);    
% QTb(n) = integral3(fun,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A),'AbsTol', 1e-6,'RelTol',1e-6);    

%Q1(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A,phi_mid(n)));    
%Q2(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A));    
%Q3(n) = integral3(fun_t,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A));    

%end
if isequal(L(1),L(2)) && isequal(M(1),M(2))
    Y1 = (r^(2*L(1) +3)./(2*L(1)+3));
else
    Y1 = 0;
end
% the differences are all numeric

Y = integral3(fun,phi_min,phi_max,theta_min,theta_max,0,r,'AbsTol', 1e-25,'RelTol',1e-25);   

Y= integral3(fun,phi_min,phi_max,theta_min,theta_max,0,r,'AbsTol', 1e-6,'RelTol',1e-6);    

