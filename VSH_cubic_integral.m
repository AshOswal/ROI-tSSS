
function Y = VSH_cubic_integral(L,M,A)

% L is the L value of the matrix
% M is the M value of the matrix
% A is the side of the cubic ROI - in metres
atol = 1e-15;
rtol = 1e-10;

if isempty(L)   
    L  = [8 8];
end
if isempty(M)
    M  = [5 5];
end
if isempty(A)
    A  = [0.1 0.9];
    %A  = 0.1;
end

% if isempty(location)
%    location = 'inner';
% end

% CONVENTION NOTE. Phi is azimuth angle measured in radians. Theta is
% elevation from negative Z axis. 

%%%%%%%%%%%%%%% TEST FOR ORTHONORMALITY OF Spherical Harmonic Functions %%%%%%%%%%%%%%

% spharm_vec.m is vectorised wrapper for spharm.m

%fun = @(theta,phi) spharm_vec(theta,phi,L(1),M(1)).*conj(spharm_vec(theta,phi,L(2),M(2))).*sin(theta);
%Q = integral2(fun,0,pi,0,2*pi,'Method','iterated','AbsTol',1e-10,'RelTol',1e-10);
%Q yields unity as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% NOTE THIS RELATES TO Vlm see appendix A1 of Taulu 2005
%fun   = @(phi,theta,r) sum(vsh_modified_in_AO(theta,phi,L(1),M(1)).*conj(vsh_modified_in_AO(theta,phi,L(2),M(2))),3) .*sin(theta).*r.^(sum(L)+2);

% actually what you want is Xlm which is given by appendix A3

fun   = @(phi,theta,r) sum(vspharm_AO(theta,phi,L(1),M(1)).*conj(vspharm_AO(theta,phi,L(2),M(2))),3) .*sin(theta).*r.^(sum(L)+2);

%fun_t = @(phi,theta,r) r.^(sum(L+2));

%% split into pyramid and upper and lower portions to make rectangle

phi_min   = [-pi/4 pi/4   3/4*pi 5/4*pi];
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


%%
for n = 1:numel(phi_min)

% the inner (deep part) from 0 -> r 
QTm(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A(1),phi_mid(n)),'AbsTol',atol,'RelTol',rtol);    
QTt(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A(1)),'AbsTol',atol,'RelTol',rtol);    
QTb(n) = integral3(fun,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A(1)),'AbsTol',atol,'RelTol',rtol);    

% the superficial part from  r -> R
QOm(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),@ (phi,theta) r(phi,theta,A(1),phi_mid(n)),A(2),'AbsTol',atol,'RelTol',rtol);    
QOt(n) = integral3(fun,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,@ (phi,theta) r_t(phi,theta,A(1)),A(2),'AbsTol',atol,'RelTol',rtol);    
QOb(n) = integral3(fun,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),@ (phi,theta) r_b(phi,theta,A(1)),A(2),'AbsTol',atol,'RelTol',rtol);    


% the inner (deep part) from 0 -> r 
% QTm1(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),0,@ (phi,theta) r(phi,theta,A(1),phi_mid(n)));    
% QTt1(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,0,@ (phi,theta) r_t(phi,theta,A(1)));    
% QTb1(n) = integral3(fun_t,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),0,@ (phi,theta) r_b(phi,theta,A(1)));    

% the superficial part from  r -> R
% QOm1(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)),@ (phi,theta) r(phi,theta,A(1),phi_mid(n)),A(2));    
% QOt1(n) = integral3(fun_t,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t,@ (phi,theta) r_t(phi,theta,A(1)),A(2));    
% QOb1(n) = integral3(fun_t,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)),@ (phi,theta) r_b(phi,theta,A(1)),A(2));        
    
end
%%
Yd      = sum(QTm+QTt+QTb);
Ys      = sum(QOm+QOt+QOb);

vs_vd      = ((4/3*pi*(A(2)^3)) - (A(1)^3))  /(A(1)^3);

Y          = vs_vd*(Yd/Ys);

% if Yd<eps||Ys <eps % a little hack 
%    Y=0;
% end

%{
if strcmp(location,'outer')
   %Y = 1./Y;
   %Y = 1-Y; % 27/9/22 change if needed back
   Y = Y; % actually change this later
   Y(isnan(Y)) = 0;
end
%}
%%
% This is essentially a sanity check showing that othogonality of VSH is obtained regardless of integrating
% over the angular dimensions of a cube 

% fun2  = @(phi,theta) sum(vsh_modified_in_AO(theta,phi,L(1),M(1)).*conj(vsh_modified_in_AO(theta,phi,L(2),M(2))),3) .*sin(theta);
% fun2  = @(phi,theta) sum(vspharm_AO(theta,phi,L(1),M(1)).*conj(vspharm_AO(theta,phi,L(2),M(2))),3) .*sin(theta);
% 
% % % % 
% for n = 1:numel(phi_min)
% Q1(n) = integral2(fun2,phi_min(n),phi_max(n),@(phi) theta_min(phi,phi_mid(n)),@(phi) theta_max(phi,phi_mid(n)));    
% Q2(n) = integral2(fun2,phi_min(n),phi_max(n),@(phi) theta_min_t(phi,phi_mid(n)),theta_max_t);    
% Q3(n) = integral2(fun2,phi_min(n),phi_max(n),theta_min_b,@(phi) theta_max_b(phi,phi_mid(n)));    
% end

