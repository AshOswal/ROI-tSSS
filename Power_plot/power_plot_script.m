clear all;
S.ROI = 'sphere';
Flmr_sum = flm_dipole_script(S);
Lin   = 11;
%
% Reconstruct Fig. 1 of the 2005 JAP paper
%
for j = 1:size(Flmr_sum{1},1)
    ff = squeeze(abs(Flmr_sum{1}(j,:,:)));
    for k = 1:size(ff,1)
        for l = 1:Lin
            indices = (l^2):((l+1)^2-1);
            ffl(k,l) = sqrt(sum(ff(k,indices).^2));
        end
    end
end
%%
A = cumsum(ffl,2);

P = A./repmat(A(:,end),1,size(A,2)); 
figure;
subplot(1,2,2);plot(P','LineWidth',2)

l   = {'0.02','0.04','0.06','0.08','0.10'}; 

xlabel('L');
ylabel('Cumulative signal power');
xlim([1 11]);ylim([0.4 1]);
set(gca,'FontSize',12);
leg = legend(l);
set(leg,'box','off');

subplot(1,2,1);plot(A','LineWidth',2);
xlabel('L');
ylabel('Signal power');
xlim([1 11]);%ylim([0 1e-2]);
set(gca,'FontSize',12);
% leg = legend(l);
% set(leg,'box','off');