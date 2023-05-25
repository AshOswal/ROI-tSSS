function [sig_out,topographies] = roi_tsss(dataset,leg)
spm('Defaults','eeg')
%MNI_roi    = [40 -20 50];%
%MNI_roi    = [0 -56 58]; % mni co-ordinate to act as centre of harmonic expansion
MNI_roi    = [0 -46 42];%
%MNI_roi    = [0 -50 40];
prefix = 'D:\home\Data\DBS-MEG\Phantom\090715\SPMstim';
inversion = 'bf';

if strcmp(dataset,'zero')
    fname = 'spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_02.mat';
elseif strcmp(dataset,'five')
    fname = 'spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_11.mat';
elseif strcmp(dataset,'twenty')
    fname = 'spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_12.mat';
elseif strcmp(dataset,'onethirty')
    fname =  'spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_07.mat';
end

%{
try
    Dn   = spm_eeg_load(fullfile(prefix,['ef' fname]));
    Dn2  = spm_eeg_load(fullfile(prefix,['efsss_' fname]));
    Dn3  = spm_eeg_load(fullfile(prefix,['efMsphere_sss_' fname]));
    Dn4  = spm_eeg_load(fullfile(prefix,['efMcube_sss_' fname]));

    % load noise covariances
    Dnoise     = spm_eeg_load('C:\home\Data\4AshO\ROI_noise.mat');
    D2noise    = spm_eeg_load('C:\home\Data\4AshO\sss_ROI_noise.mat');
    weight     = load(fullfile(prefix,[dataset '_weight.mat']));
    weight     = weight.weight;
    weight1    = load(fullfile(prefix,[dataset '_weight1.mat']));
    weight1    = weight1.weight1;

    
catch
%}
    D  = spm_eeg_load([ 'D:\home\Data\DBS-MEG\Phantom\090715\SPMstim\' fname]);
    %D          = spm_eeg_load('V:\vlad_shared\4AshO\ROI_sim.mat');
    %Dnoise     = spm_eeg_load('C:\home\Data\4AshO\ROI_noise.mat');
    
    D = forward_model_phantom(D);
    
    %% convert from MNI to headspace then to the TSSS space
    
    trans           = D.inv{1}.forward.fromMNI;
    head_roi_ctf    = spm_eeg_inv_transform_points(trans,MNI_roi);  % CTF co-ordinates
    fid             = D.fiducials;
    fid             = ft_convert_units(fid, 'm');
    
    %Mc = spm_eeg_inv_headcoordinates(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
    
    M = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
    head_roi_elekta = spm_eeg_inv_transform_points(M,head_roi_ctf); % Centre of harmonic expansion
    
    %% apply raw first
    S            = [];
    S.D          = D;
    S.tsss       = 1;
    S.t_window   = 5;
    S.xspace     = 0;
    D2           = tsss_spm_enm(S);
    %{
    S.D          = Dnoise;
    D2noise      = tsss_spm_enm(S);
    %}
    %% process the original data through a new ROI-tSSS pipeline for spherical or cubic regions of interest
    S                = [];
    S.D              = D2;
    S.r_sphere       = head_roi_elekta';
    S.Lin            = 8;
    S.roi_rad        = 0.03;
    S.corr_limit     = 1;
    S.input          = 'inner';
    S.roi_type       = 'sphere';
    S.prefix         = 'ROI_tsss_sphere_';
    S.interference_space = 'extern_';
    S.temporal_samples_threshold = 2.5e4;
    S.apply_temporal = 1;
    % new option to normalise weights 
    S.normalise      = 1;
    [D3,~,weight]      = tsss_roi_ext(S);
    save(fullfile(prefix,[dataset '_' 'weight.mat']),'weight');
    D3 = forward_model_phantom(D3);
    %{
    S.D              = D2noise;
    S.prefix         = 'ROI_tsss_sphere_noise_';
    [D3noise,~]      = tsss_roi_ext(S);
    %}
    %%
   
    %{
    S.D              = D2;
    S.roi_type       = 'cube';
    S.prefix         = 'ROI_tsss_cube_';
    %S.apply_temporal = 0;
    %S.corr_limit     = 1;
    [D4,~,weight1]     = tsss_roi_ext(S);
    save(fullfile(prefix,[dataset '_' 'weight1.mat']),'weight1');
    D4 = forward_model_phantom(D4);
    %}
    % the noise dataset
    %{
    S.D              = D2noise;
    S.prefix         = 'ROI_tsss_cube_noise';
    [D4noise,~]      = tsss_roi_ext(S);
    %}
    %%
    S = [];
    S.D = D3;
    S.task = 'project3D';
    S.modality = 'MEG';
    S.updatehistory = 0;
    S.save    = 1;
    D3 = spm_eeg_prep(S);
    %S.D = D4;
    %D4 = spm_eeg_prep(S);
   
    S = [];
    S.band = 'high';
    S.freq = 5;
    S.D = D;
    D = spm_eeg_filter(S);
    S.D = D2;
    D2 = spm_eeg_filter(S);
    S.D = D3;
    D3 = spm_eeg_filter(S);
    %S.D = D4;
    %D4 = spm_eeg_filter(S);
   
        
    nsamples      = floor(4*D3.fsample/100); % length of trials in samples for epoching
    
    % deal with the simulated LFP signal and partition this according
    sim_dip= ft_preproc_bandpassfilter(squeeze(D3(D3.indchannel('LFP_CLEAN'),:,1)),D3.fsample,[26 28]);

    sim_dip_smooth = smooth(sim_dip,10);
    
    [pks,locs] = findpeaks(sim_dip_smooth,'MinPeakDistance',floor(D3.fsample/27)-1);
    ns         = min(unique(diff(locs)));
        
    %trl    = [locs locs+nsamples zeros(size(locs))];
    trl   = [locs locs+ns-1 zeros(size(locs))];    
    out   = find(trl(:,2)>numel(sim_dip_smooth));     
    trl(out,:) = [];
    locs(out) = [];  
    
    i = 1;
    while i<numel(locs) && any(locs(i+1:end)-locs(i) <= nsamples)
        %del = find(locs(i+1:end)-locs(i) <= nsamples);
        del = find(locs(i+1:end)-locs(i) <= ns);
        
        trl(del+i,:) = [];
        locs(del+i)  = [];
        i = i+1;
    end
    average_evoked = [];
    for n = 1:size(trl,1)
        average_evoked = [average_evoked;sim_dip_smooth(trl(n,1):trl(n,2))'];
    end
    average_evoked = mean(average_evoked,1);
    average_evoked = (average_evoked-mean(average_evoked))./std(average_evoked);
    
    
    S    = [];
    S.D  = D;
    S.trl = trl;
    S.bc = 0;
    Dn = spm_eeg_epochs(S);
    
    S.D  = D2;
    Dn2 = spm_eeg_epochs(S);
    
    S.D  = D3;
    Dn3 = spm_eeg_epochs(S);
   
    %S.D  = D4;
    %Dn4 = spm_eeg_epochs(S);
    
%end

%% noise covariances

%labels = Dn.chanlabels(Dn.indchantype('MEG'));
%ic     = indchannel(Dnoise,labels);


% Dnoise_cov  = cov(Dnoise(:,:)');
% D2noise_cov = cov(D2noise(:,:)');
% D3noise_cov = cov(D3noise(:,:)');%weight{1}*weight{1}';
% D4noise_cov = cov(D4noise(:,:)');%weight1{1}*weight1{1}';


Dnf   = ftraw(Dn,setdiff(Dn.indchantype('MEG'),Dn.badchannels),1:50,':');
Dn2f  = ftraw(Dn2,setdiff(Dn2.indchantype('MEG'),Dn2.badchannels),1:50,':');
Dn3f  = ftraw(Dn3,setdiff(Dn3.indchantype('MEG'),Dn3.badchannels),1:50,':');
%Dn4f  = ftraw(Dn4,setdiff(Dn4.indchantype('MEG'),Dn4.badchannels),1:50,':');

Dna   = ftraw(Dn,setdiff(Dn.indchantype('MEG'),Dn.badchannels),':',':');
Dn2a  = ftraw(Dn2,setdiff(Dn2.indchantype('MEG'),Dn2.badchannels),':',':');
Dn3a  = ftraw(Dn3,setdiff(Dn3.indchantype('MEG'),Dn3.badchannels),':',':');
%Dn4a  = ftraw(Dn4,setdiff(Dn4.indchantype('MEG'),Dn4.badchannels),':',':');

cfg            = [];
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'no';
Dnf             = ft_timelockanalysis(cfg,Dnf);
Dn2f            = ft_timelockanalysis(cfg,Dn2f);
Dn3f            = ft_timelockanalysis(cfg,Dn3f);
%Dn4f            = ft_timelockanalysis(cfg,Dn4f);
Dna             = ft_timelockanalysis(cfg,Dna);
Dn2a            = ft_timelockanalysis(cfg,Dn2a);
Dn3a            = ft_timelockanalysis(cfg,Dn3a);
%Dn4a            = ft_timelockanalysis(cfg,Dn4a);
%%
cfg             = [];
cfg.parameter   = 'avg';
cfg.colorbar    = 'eastoutside';
cfg.style      = 'straight';
cfg.marker      = 'off';
%cfg.interactive = 'no';
cfg.figure      = 'no';
% if strcmp(dataset,'zero')
% cfg.zlim       = [-25 25];
% end
cfg.comment    = 'no';
%cfg.colorbar   = 'southoutside';

figure;
ft_topoplotER(cfg,Dnf);
title('Standard');
set(gca,'FontSize',22);

figure;
ft_topoplotER(cfg,Dn2f)
title('tSSS');
set(gca,'FontSize',22);

figure;
ft_topoplotER(cfg,Dn3f)
set(gca,'FontSize',22);
title('ROI-tSSS sphere')

% figure;
% ft_topoplotER(cfg,Dn4f)
% set(gca,'FontSize',22);
% title('ROI-tSSS cube')

%%

switch inversion
    case 'dipole'
        data    = spm_eeg_inv_get_vol_sens(Dn,1,'MNI-aligned','inv','MEG');
        grad    = data.MEG.sens;

        data2   = spm_eeg_inv_get_vol_sens(Dn2,1,'MNI-aligned','inv','MEG');
        grad2   = data2.MEG.sens;

        data3   = spm_eeg_inv_get_vol_sens(Dn3,1,'MNI-aligned','inv','MEG');
        grad3   = data3.MEG.sens;

        %data4   = spm_eeg_inv_get_vol_sens(Dn4,1,'MNI-aligned','inv','MEG');
        %grad4   = data4.MEG.sens;

        cfg=[];
        cfg.appenddim = 'chan';
        cfg.headmodel = data.MEG.vol;
        cfg.grad      = grad;
        cfg.grid.resolution = 0.01;
        cfg.gridsearch = 'yes';
        cfg.channel = 'all';
        cfg.model  = 'regional';
        %cfg.reducerank = 'no';
        %cfg.dipfit.noisecov = Dnoise_cov;
        %cfg.weight = [];

        % standard
        cfg.headmodel = data.MEG.vol;
        cfg.grad      = grad;
        dip_source   = ft_dipolefitting(cfg, Dna);
        % tSSS
        cfg.headmodel = data2.MEG.vol;
        cfg.grad      = grad2;
        %cfg.dipfit.mleweight = inv(D2noise_cov);
        dip_source2  = ft_dipolefitting(cfg, Dn2a);
        % tsss + ROIs
        cfg.headmodel = data3.MEG.vol;
        cfg.grad      = grad3;
        % cfg.dipfit.noisecov = D3noise_cov;
        % cfg.dipfit.mleweight = pinv(D3noise_cov);
        dip_source3  = ft_dipolefitting(cfg, Dn3a);
        % tsss + ROIc
        %cfg.headmodel = data4.MEG.vol;
        %cfg.grad      = grad4;
        % cfg.dipfit.noisecov = D4noise_cov;
        %dip_source4  = ft_dipolefitting(cfg, Dn4a);

        %% source estimation

        [u, s, v] = svd(dip_source.dip.mom);
        est_ori = u(:,1)'; % orientation
        signal = s(1,1)*v(:,1)';

        [u, s, v] = svd(dip_source2.dip.mom);
        est_ori2 = u(:,1)'; % orientation
        signal2 = s(1,1)*v(:,1)';

        [u, s, v] = svd(dip_source3.dip.mom);
        est_ori3 = u(:,1)'; % orientation
        signal3 = s(1,1)*v(:,1)';

        %[u, s, v] = svd(dip_source4.dip.mom);
        %est_ori4 = u(:,1)'; % orientation
        %signal4 = s(1,1)*v(:,1)';

        %signal5 = fliplr(signal5)est_ori*dip_source5.dip.mom;%

      
    case 'bf'


        matlabbatch = {};
        matlabbatch{1}.spm.tools.beamforming.data.dir = {'D:\home\Data\4AshO\test'};
        matlabbatch{1}.spm.tools.beamforming.data.val = 1;
        matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
        matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
        matlabbatch{1}.spm.tools.beamforming.data.overwrite = 1;
        matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
        matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
        matlabbatch{2}.spm.tools.beamforming.sources.plugin.mni_coords.pos = MNI_roi;
        matlabbatch{2}.spm.tools.beamforming.sources.normalise_lf = false;
        matlabbatch{2}.spm.tools.beamforming.sources.visualise = true;
        matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
        matlabbatch{3}.spm.tools.beamforming.features.woi = [-Inf Inf];
        matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEG'};
        matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
        matlabbatch{3}.spm.tools.beamforming.features.cross_terms = 'megeeg';
        matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.foi = [0 45];
        matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.taper = 'hanning';
        matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 5;
        matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
        matlabbatch{3}.spm.tools.beamforming.features.visualise = 0;
        matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
        matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{5}.spm.tools.beamforming.output.plugin.montage.method = 'max';
        matlabbatch{5}.spm.tools.beamforming.output.plugin.montage.vois = cell(1, 0);
        matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));

        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'write';
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'MEG';
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.none = 0;
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.prefix = 'B';

        %{
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'online';
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'MEG';
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.none = 0;
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.prefix = 'B';
        %}
        matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(Dn)};
        spm_jobman('run',matlabbatch);

        matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(Dn2)};
        spm_jobman('run',matlabbatch);

        matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(Dn3)};
        spm_jobman('run',matlabbatch);

        %matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(Dn4)};
        %spm_jobman('run',matlabbatch);

        [a,~] = fileparts(fullfile(Dn));
        dn    = spm_eeg_load(fullfile(a,['B' Dn.fname]));
        dn2   = spm_eeg_load(fullfile(a,['B' Dn2.fname]));
        dn3   = spm_eeg_load(fullfile(a,['B' Dn3.fname]));
        %dn4   = spm_eeg_load(fullfile(a,['B' Dn4.fname]));

        signal = mean(dn(:,:,:),3);
        signal2 = mean(dn2(:,:,:),3);
        signal3 = mean(dn3(:,:,:),3);
        %signal4 = mean(dn4(:,:,:),3);


end

figure;
plot(Dn3.time,(signal-mean(signal))./std(signal),'b','LineWidth',2);hold on;
plot(Dn3.time,(signal2-mean(signal2))./std(signal2),'r--','LineWidth',2);hold on;
plot(Dn3.time,(signal3-mean(signal3))./std(signal3),'m--','LineWidth',2);hold on;
%plot(Dn3.time,(signal4-mean(signal4))./std(signal4),'g--','LineWidth',2);hold on;

xlim([0 0.035]);
ylim([-2 2]);
set(gca,'FontSize',20);
xlabel('Time (s)');
ylabel('Source Amplitude (AU)');

xlim([0 0.035]);
%ylim([-2 2]);
set(gca,'FontSize',20);
xlabel('Time (s)');
ylabel('Source Amplitude (AU)');

if leg
    h = legend('Standard','tSSS','ROI-tSSS sphere','ROI-tSSS cube');
    set(h,'box','off','location','EastOutside');
end

sig_out      = [signal;signal2;signal3;];%signal4];
topographies =  struct('nil',Dnf.avg,'tsss',Dn2f.avg,'sroi_tsss',Dn3f.avg);%,'croi_tsss',Dn4f.avg);
%%

function D = forward_model_phantom(D)

% Note this bit is only for the phantom
% ---- TO DO - add code to get the forward model from the data and reapply
% this to the new montage. 

clear batch
job1{1}.spm.meeg.source.headmodel.D = {fullfile(D)};
job1{1}.spm.meeg.source.headmodel.val = 1;
job1{1}.spm.meeg.source.headmodel.comment = '';
job1{1}.spm.meeg.source.headmodel.meshing.meshres = 2;

job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'FIL_CTF_L';
job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'FIL_CTF_R';

job1{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
job1{1}.spm.meeg.source.headmodel.forward.meg = 'Single Sphere';

spm_jobman('run', job1);

D = reload(D);

% Correction of the head model

M = eye(4);
M(1:3, 4) = -1e3*D.inv{1}.forward.vol.o';

D.inv{1}.datareg(1).fid_mri = ft_transform_headshape(M, D.inv{1}.datareg(1).fid_mri);
D.inv{1}.datareg(1).toMNI = D.inv{1}.datareg(1).toMNI/M;
D.inv{1}.datareg(1).fromMNI = inv(D.inv{1}.datareg(1).toMNI);
D.inv{1}.datareg(1).modality = 'MEG';

spm_eeg_inv_checkdatareg(D);


D = spm_eeg_inv_forward(D);

D.inv{1}.forward.vol.o = [0 0 0];
D.inv{1}.forward.vol.r = 0.075;

spm_eeg_inv_checkforward(D);

save(D);


%{
[vol, sens] = ft_prepare_vol_sens(data.MEG.vol, data.MEG.sens, 'channel', D.chanlabels(D.indchantype('MEG')));
lf  = ft_compute_leadfield(dip_source.dip.pos, sens, vol, 'dipoleunit', 'nA*m', 'chanunit', D.units(D.indchantype('MEG')));

[vol, sens] = ft_prepare_vol_sens(data1.MEG.vol, data1.MEG.sens, 'channel', D3.chanlabels(D3.indchantype('MEG')));
lf2 = ft_compute_leadfield(dip_source3.dip.pos, sens, vol, 'dipoleunit', 'nA*m', 'chanunit', D.units(D.indchantype('MEG')));

lft = lf*est_ori';
lft2= lf2*est_ori';

Dn_s.var(:,2) = [];
Dn_s.time(2)  = [];
Dn_s.dof(:,2) = [];
a = Dn_s;
b = Dn_s;
a.avg = lft;
b.avg = lft2;

cfg            = [];
cfg.parameter  = 'avg';
cfg.colorbar   = 'east';
cfg.style      = 'straight';
cfg.marker     = 'off';
%cfg.zlim       = [-5 5];
cfg.comment    = 'no';

figure;
subplot(2,1,1)
title('unprocessed');
ft_topoplotTFR(cfg,a);
subplot(2,1,2)
title('SSS');
ft_topoplotTFR(cfg,b)

%%

figure;
%subplot(1,2,1)
plot(Dn3.time,(signal-mean(signal))./std(signal),'r','LineWidth',2);hold on;
plot(Dn3.time,(signal3-mean(signal3))./std(signal3),'k','LineWidth',2);hold on;
set(gca,'FontSize',14);
xlabel('Time (s)');
ylabel('Source Amplitude (AU)');
legend('Standard pre-processing','Standard pre-processing + SSS-ROI');
legend boxoff

%%
plot([(signal-mean(signal))./std(signal);(signal2-mean(signal2))./std(signal2);(signal3-mean(signal3))./std(signal3);average_evoked]')
%axis([0 500 -1.5e7 1.5e7])
legend('unprocessed','SSS','SSS+ROI')
%load('est_ori');load('est_pos');
%%
figure
ft_plot_dipole(dip_source.dip.pos, est_ori, 'color', 'b' );hold on
ft_plot_dipole(dip_source2.dip.pos, est_ori2, 'color', 'g' );hold on
ft_plot_dipole(dip_source3.dip.pos, est_ori3, 'color', 'r' );hold on
%ft_plot_dipole(est_pos,est_ori,'color','k');
rotate3d off;
axis off
axis vis3d


figure;
plot([(signal3-mean(signal3))./std(signal3);average_evoked]')
%axis([0 500 -1.5e7 1.5e7])

%}