function [G,sup_scale,deep_scale] = GsGd(L,N,R,Norm)

% this is the matrix which will be multiplied by the alpha coefficients in
% order to get data for the superficial or deep components.

if isempty(L)
   L = 7;
end
if isempty(N)
   N = 7;
end
if isempty(R)
   R = [0.06 0.17];
end
if isempty(Norm)
   Norm = 0;
end

%%

% we must also normalise the filter weight. Note by construction that this
% will be 0 when r -> 0. We will also introduce an approximation in that
% the filter weight should be 1 if the square fits inside the sphere; this
% derives using Pythagoras' theorem

%side = ((4/3)*pi*(R(2)^3))^(1/3);
Rn = [(sqrt(3)/2)*R(2) R(2)];
%%
G = [];
Gn= [];
tic

L1s = [];
for L1 = 1:L
    parfor L2 = -L1:L1 
        y = compute_integral_G([L1 L2],N,R);
        %y = y*((2*L1)+3)/3;
        % novel eigenvalue scaling
        L1s = [L1s,L1];
        %y = y*(2*L1+3)./5;
        G = [G;y];
        
        % this was the normaliser but actually we can just take Tolga's
        %yn= compute_integral_G([L1 L2],N,Rn,T);
        %Gn= [Gn;yn];
    end
end
toc

%%

if Norm
    deep_scale  = diag((2*L1s+3)./5);
    sup_scale   = diag((2*L + 3)./(2*L1s+3));
else
    deep_scale  = diag((2*L1s+3)./3);
    sup_scale   = [];
end

%%    

%G = G/Gn;
%G = G*(L(1)+L(2))/3;
if 0
figure;
subplot(1,3,1);imagesc(G);title('r');
subplot(1,3,2);imagesc(Gn);title('R');
subplot(1,3,3);imagesc(G/Gn);title('Normalised');
end
end


function y = compute_integral_G(L,N,R)
% semi-recursive routine

y = [];
for n1=-N:N
    % symmetric matrix for expediency just take the upper triangle/diagonal
    if ~isequal(L,[N n1]) ||  N<L(1) || (N==L(1) && n1<L(2))     % last bit takes diagonal - it is sufficient alone
        y = [y 0];
    else
        y = [y VSH_cubic_integral([L(1) N],[L(2),n1],R)];
    end
end
% recursive call
if N>1
    y = [compute_integral_G(L,N-1,R)  y];
end


end