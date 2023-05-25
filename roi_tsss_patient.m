clear all
eeglab;
%spmfile = 'D:\home\Data\DBS-MEG\LN_C21\170314\SPMstim\LN_C21_STIM_R_20_off.mat';
spmfile     = 'D:\home\Data\DBS-MEG\LN_C10\270613\SPMstim\LN_C10_STIM_R_5_off.mat';
stimfrq     = 5;
peak_offset = 50;

D       = spm_eeg_load(spmfile);

%% apply tsss

S            = [];
S.D          = D;
S.tsss       = 1;
S.corr_limit = 0.95;
S.t_window   = 5;
S.xspace     = 0;
D_tsss       = tsss_spm_enm(S);

%%


MNI_roi         = [40 -20 50];
trans           = D.inv{1}.forward.fromMNI;
head_roi_ctf    = spm_eeg_inv_transform_points(trans,MNI_roi);  % CTF co-ordinates
fid             = D.fiducials;
fid             = ft_convert_units(fid, 'm');

%Mc = spm_eeg_inv_headcoordinates(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));

M = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
head_roi_elekta = spm_eeg_inv_transform_points(M,head_roi_ctf); % Centre of harmonic expansion

S                = [];
S.D              = D_tsss;
S.r_sphere       = head_roi_elekta';
S.Lin            = 8;
S.roi_rad        = 0.03;
S.corr_limit     = 0.98;
S.input          = 'inner';
S.roi_type       = 'sphere';
S.prefix         = 'ROI_tsss_sphere_';
S.interference_space = 'extern_';
S.temporal_samples_threshold = 2.5e4;
S.apply_temporal = 1;
% new option to normalise weights
S.normalise      = 1;
[D,~,weight]     = tsss_roi_ext(S);
% save(fullfile(prefix,[dataset '_' 'weight.mat']),'weight');
% D = forward_model_phantom(D);


%%

nsamples      = floor(4*D.fsample/100); % length of trials in samples for epoching

%[t1, t2, ind] = stim_bounds(squeeze(D(D.indchannel('STIM'),:,1)), D.fsample, 5)

% deal with the simulated LFP signal and partition this according
sim_dip = squeeze(D(D.indchannel('STIM'),:,1))';


[pks,locs] = findpeaks(diff(sim_dip),'MinPeakDistance',floor(D.fsample/stimfrq)-3);
for i = 1:numel(locs)
    ind = 0;
    while sim_dip(locs(i)-ind) > sim_dip(locs(i) - (ind+1))
        ind = ind+1;
    end
    pre(i) = ind;
end

ns         = min(unique(diff(locs)));
locs       = locs-pre';
trl        = [locs locs+ns-1 zeros(size(locs))];
out        = [find(trl(:,2)>numel(sim_dip))  find(trl(:,1)<5)];
trl(out,:) = [];
locs(out)  = [];

i = 1;
while i<numel(locs) && any(locs(i+1:end)-locs(i) <= nsamples)
    del = find(locs(i+1:end)-locs(i) <= ns);
    
    trl(del+i,:) = [];
    locs(del+i)  = [];
    i = i+1;
end

%% trl adjustment
trl(:,[1 2]) = trl(:,[1 2])-round(D.fsample*peak_offset/1e3);

%%
average_evoked = [];
for n = 1:size(trl,1)
    average_evoked = [average_evoked;sim_dip(trl(n,1):trl(n,2))'];
end
average_evoked = mean(average_evoked,1);
average_evoked = (average_evoked-mean(average_evoked))./std(average_evoked);


S           = [];
S.D         = D;
S.trl       = trl;
S.bc        = 0;
Dn          = spm_eeg_epochs(S);
Dn          = timeonset(Dn,-peak_offset*1e-3);
save(Dn);
S.D         = D_tsss;
Dn_tsss     = spm_eeg_epochs(S);
Dn_tsss     = timeonset(Dn_tsss,-peak_offset*1e-3);

%% baseline correct

S           = [];
S.D         = Dn;
S.timewin   = [-50 -30];
%S.timewin   = [-20 -5];
Dn          = spm_eeg_bc(S);
S.D         = Dn_tsss;
Dn_tsss     = spm_eeg_bc(S);

%%

S          = [];
S.D        = Dn;
S.channels = D.chanlabels(D.indchantype('MEG'));
%S.timewin = [-0.02 0.05]*1e3;
S.timewin  = [-0.02 0.02]*1e3;
mChan      = spm_eeg_crop(S);
S.D        = Dn_tsss;
mChan_tsss = spm_eeg_crop(S);

mChan_sd      = std(squeeze(mean(mChan(:,:,:),1)),1,2)./sqrt(mChan.nchannels);
mChan_tsss_sd = std(squeeze(mean(mChan_tsss(:,:,:),1)),1,2)./sqrt(mChan_tsss.nchannels);

% add a montage
mon          = [];
mon.tra      = ones(1,size(mChan,1))./size(mChan,1);
mon.labelnew = {'average_ch'};
mon.labelorg = mChan.chanlabels;
mChan        = montage(mChan,'add',mon);
mChan        = montage(mChan,'switch',0);

S             = [];
S.D           = mChan;
mChan_av      = spm_eeg_average(S);
S.D           = mChan_tsss;
mChan_tsss_av = spm_eeg_average(S);

%% SPM images
% S      = [];
% S.D    = mChan_av;
% S.mode = 'scalp x time';
% %S.timewin = [3 200]*1e-3;
% iname          = spm_eeg_convert2images(S);
% 
% [pre,post] = spm_fileparts(iname);
% imname     = fullfile(pre,[post '.nii']);
% sname      = fullfile(pre,['s' [post '.nii']]);
% spm_smooth(imname,sname,[2.5 2.5 0.5])
%%
amp = mean(mChan_av(:,:),1);
pk  = spm_percentile(amp,98); 
[pks,locs] = findpeaks(amp,'MinPeakHeight',pk);

figure;
plot(mChan_av.time,mean(mChan_av(:,:),1),'k','LineWidth',2);hold on;
plot(mChan_av.time(locs),pks,'bs','MarkerSize',8,'MarkerFaceColor','b')
for l=1:numel(locs)
    str = [sprintf('%.1f',mChan_av.time(locs(l))*1e3) 'ms'];
    text(mChan_av.time(locs(l)),pks(l)+4,str,'FontSize',14,'Color','b')
end
set(gca,'FontSize',14);hold on
% add standard eror
fill([mChan_av.time';flipud(mChan_av.time')],[mean(mChan_av(:,:),1)'-mChan_sd;flipud(mean(mChan_av(:,:),1)'+mChan_sd)],[0.2 0.2 0.2],'linestyle','none','FaceAlpha',0.3);
% hold on 
% plot(mChan_av.time,mChan_sd) n ku/, 
yo = ylim; 
%ylim([yo(1) yo(2)+10]);
xlabel('Time (s)');ylabel('Relative Amplitude');

figure;
plot(mChan_tsss_av.time,mean(mChan_tsss_av(:,:),1),'k','LineWidth',2);hold on;

%%
erp = squeeze(mean(Dn(Dn.indchantype('MEG'),:,:),1));
figure;
[outdata, outvar, outtrl] = erpimage(erp,[],Dn.time*1e3,[],15,[],'avg_type','Gaussian','phasesort',[2.4 0 1/diff(mChan.time(locs))],'cycles',4);
figure;
[outdata2, outvar2, outtrl2] = erpimage(erp,[],Dn.time*1e3,[],15,[],'avg_type','Gaussian');
%%
figure;
timeind = abs(Dn.time)<=0.01;
subplot(2,1,1);
imagesc(Dn.time(timeind)*1e3,outtrl, outdata(timeind,:)',[-75 75])
set(gca,'YDir','normal','FontSize',14);
xlabel('Time ms');ylabel('Sorted Trials');
hcb=colorbar;
hcb.Title.String = "Amplitude";

subplot(2,1,2);
imagesc(Dn.time(timeind)*1e3,outtrl2, outdata2(timeind,:)',[-75 75])
set(gca,'YDir','normal','FontSize',14);
xlabel('Time ms');ylabel('Trials')
hcb=colorbar;
hcb.Title.String = "Amplitude";
%%

for l=1:numel(locs)
    Dnf            = ftraw(mChan_av,setdiff(mChan_av.indchantype('MEG'),mChan_av.badchannels),[locs(l)-1 locs(l) locs(l)+1],':');
    Dnf_tsss       = ftraw(mChan_tsss_av,setdiff(mChan_tsss_av.indchantype('MEG'),mChan_tsss_av.badchannels),[locs(l)-1 locs(l) locs(l)+1],':');
    cfg            = [];
    cfg.channel    = 'all';
    cfg.trials     = 'all';
    cfg.keeptrials = 'no';
    Dnf            = ft_timelockanalysis(cfg,Dnf);
    Dnf_tsss       = ft_timelockanalysis(cfg,Dnf_tsss);



    cfg             = [];
    cfg.parameter   = 'avg';
    %cfg.colorbar    = 'eastoutside';
    cfg.style      = 'straight';
    cfg.marker      = 'off';

    %cfg.interactive = 'no';
    cfg.figure      = 'no';
    % if strcmp(dataset,'zero')
    cfg.zlim       = [-200 1500];
    % end
    cfg.comment    = 'no';
    cfg.colorbar   = 'southoutside';

    figure;
    ft_topoplotER(cfg,Dnf);
    title([sprintf('%.1f',mChan_av.time(locs(l))*1e3) ' ms']);
    set(gca,'FontSize',20);

%     figure;
%     ft_topoplotER(cfg,Dnf_tsss);
end