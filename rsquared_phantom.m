% rsquared_phantom
clear all;close all
[sig_out,topographies]   = roi_tsss('zero',1);
[sig_out2,topographies2] = roi_tsss('five',1);
[sig_out3,topographies3] = roi_tsss('twenty',1);
[sig_out4,topographies4] = roi_tsss('onethirty',1);


%%

error_stan_5      = mean(mean((topographies2.nil - topographies.nil).^2,2));
error_tsss_5      = mean(mean((topographies2.tsss - topographies.tsss).^2,2));
error_sroi_5      = mean(mean((topographies2.sroi_tsss - topographies.sroi_tsss).^2,2));
error_croi_5      = mean(mean((topographies2.croi_tsss - topographies.croi_tsss).^2,2));

error_stan_20      = mean(mean((topographies3.nil - topographies.nil).^2,2));
error_tsss_20      = mean(mean((topographies3.tsss - topographies.tsss).^2,2));
error_sroi_20      = mean(mean((topographies3.sroi_tsss - topographies.sroi_tsss).^2,2));
error_croi_20      = mean(mean((topographies3.croi_tsss - topographies.croi_tsss).^2,2));

error_stan_130      = mean(mean((topographies4.nil - topographies.nil).^2,2));
error_tsss_130      = mean(mean((topographies4.tsss - topographies.tsss).^2,2));
error_sroi_130      = mean(mean((topographies4.sroi_tsss - topographies.sroi_tsss).^2,2));
error_croi_130      = mean(mean((topographies4.croi_tsss - topographies.croi_tsss).^2,2));

error5  = [error_stan_5,error_tsss_5,error_sroi_5,error_croi_5];
error20  = [error_stan_20,error_tsss_20,error_sroi_20,error_croi_20];
error130  = [error_stan_130,error_tsss_130,error_sroi_130,error_croi_130];
%%
figure;
subplot(2,1,1);
msize = 10;
plot(1,log(error_stan_5),'rs','MarkerSize',msize,'MarkerFaceColor','r');hold on
plot(1,log(error_stan_20),'b^','MarkerSize',msize,'MarkerFaceColor','b');hold on
plot(1,log(error_stan_130),'ko','MarkerSize',msize,'MarkerFaceColor','k');
xlim([0.75 4.25]);
set(gca,'XTickLabels',[],'FontSize',14);
ylabel('Log MSE')
h = legend('5Hz','20Hz','130Hz');
set(h,'box','off','location','Northeast');

subplot(2,1,2);
plot(2,error_tsss_5,'rs','MarkerSize',msize,'MarkerFaceColor','r');hold on
plot(2,error_tsss_20,'b^','MarkerSize',msize,'MarkerFaceColor','b');hold on
plot(2,error_tsss_130,'ko','MarkerSize',msize,'MarkerFaceColor','k');hold on

plot(3,error_sroi_5,'rs','MarkerSize',msize,'MarkerFaceColor','r');hold on
plot(3,error_sroi_20,'b^','MarkerSize',msize,'MarkerFaceColor','b');hold on
plot(3,error_sroi_130,'ko','MarkerSize',msize,'MarkerFaceColor','k');hold on

plot(4,error_croi_5,'rs','MarkerSize',msize,'MarkerFaceColor','r');hold on
plot(4,error_croi_20,'b^','MarkerSize',msize,'MarkerFaceColor','b');hold on
plot(4,error_croi_130,'ko','MarkerSize',msize,'MarkerFaceColor','k');hold on

xlim([0.75 4.25]);
set(gca,'FontSize',14,'XTick',[1,2,3,4],'XTickLabel',{'Standard','tSSS','ROI-tSSS sphere','ROI-tSSS cube'});
ylabel('MSE');




