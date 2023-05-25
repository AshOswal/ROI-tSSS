function [flm,indices] = flm_dipole(r,rq,q,Lin)

% r  - r distance
% rq - position of dipole
% q  - current
rqn = norm(rq);
theta = acos(rq(3)/rqn);
phi = atan2(rq(2),rq(1));

% qs(1) - unit vector e(theta)/ qs(2) - unit vector e(phi)
qs(1) = q(1)*cos(theta)*cos(phi) + q(2)*cos(theta)*sin(phi) - q(3)*sin(theta);
qs(2) = -q(1)*sin(phi) + q(2)*cos(phi);
% only theta and phi components are needed as Xlm depends only on these

count = 1;
for l = 1:Lin
   for m = -l:l
      %disp([l m]); 
      indices(count,:) = [l m];
      flm(count) = i*sqrt(l/(2*l+1))*(rqn^(l+2))*dot(conj(vspharm(theta,phi,l,m)),qs)/r^(l+2);
      %flm(count) = (-i/(2*l+1))*sqrt(l/(l+1))*(rqn^l)*dot(qs,vspharm(theta,phi,l,m))/r^(l+2);
      
      % we need to multiply this by G for different ROI radii
      
      count = count + 1;
   end
end

