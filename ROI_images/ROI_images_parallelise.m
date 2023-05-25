% code for making a beamformer image from regions of interest

clear all;
spm('defaults','eeg');

%% subject settings
subject                          = 'LN_C10';
[files, sequence, root, details] =dbs_subjects(subject, 0);

image_name          = {fullfile(root,'SPMstim/evoked_R_5_i1.nii');...
    fullfile(root,'SPMstim/evoked_R_5_i2.nii');...
    };

%image_name          = 'C:\home\Data\4AshO\ROI_images\power_image.nii';
%spmfile             = 'C:\home\Data\4AshO\ROI_images\ROI_sim.mat';
%spmfile             = fullfile(root,'SPMstim\beROI_tsss_sphere_einner_rsss_LN_C10_STIM_R_5_off.mat');
spmfile             = fullfile(root,'SPMstim\sss_LN_C10_STIM_R_5_off.mat');
bf_dir              = 'D:\home\Data\4AshO\ROI_images';

%% epoch settings
stimfrq             = 5;
peak_offset         = 50;

%% configuration settings
forward_model_space = 'MNI-aligned';  % space of forward model
grid_space          = 'MNI template'; % space of grid
image_space         = 'MNI template';
inv_index           = 1;
bound               = 'iskull';


% beamformer settings
resolution          = 10;
visualise_sources   = 1;
timewin             = [-50 120];%[3 150];
freqoi              = [0 Inf];
frqout              = [-Inf Inf];

ROI                 = 1;
normpow             = 0;

% output settings
% - this will either be power, RMS of timeseries over timewin or peakamp at
% a peaktime(s)
output.type         = 'peakamp';
output.time         =  [2.5 4.2];

side                = 'R'; %{'R','L','B'}
%%
D = spm_eeg_load(spmfile);

if ~isfield(D,'inv')
    warning('running headmodel as this has not been specified');
    D = dbs_meg_headmodelling(D);
end

Data           = spm_eeg_inv_get_vol_sens(D,1,'MNI-aligned','inv','MEG');

if strcmp(bound,'iskull')
    constraint     = export(gifti(D.inv{inv_index}.mesh.tess_iskull), 'ft');
elseif strcmp(bound,'scalp')
    constraint     = export(gifti(D.inv{inv_index}.mesh.tess_scalp), 'ft');
end

M1 = Data.transforms.toNative;

switch grid_space
    case 'MNI template'
        M1 = Data.transforms.toMNI/M1;
        M2 = inv(Data.transforms.toMNI);
    case 'MNI-aligned'
        M1 = Data.transforms.toMNI_aligned/M1;
        M2 = inv(Data.transforms.toMNI_aligned);
    case 'Head'
        M1 = Data.transforms.toHead/M1;
        M2 = inv(Data.transforms.toHead);
    case 'Native'
        M2 = inv(M1);
        M1 = eye(4);
end

constraint = ft_determine_units(ft_transform_geometry(M1, constraint));

mn = min(constraint.pnt);
mx = max(constraint.pnt);

if isequal(constraint.unit, 'm')
    resolution = 1e-3*resolution;
end

% If zero is inside the brain, make sure grid points fall on multiples of
% resolution to ease simulating data from points on the grid
if mn(1)<0 && mx(1)>0
    grid.xgrid = [fliplr(0:-resolution:mn(1)) resolution:resolution:mx(1)];
else
    grid.xgrid = mn(1):resolution:mx(1);
end

if mn(2)<0 && mx(2)>0
    grid.ygrid = [fliplr(0:-resolution:mn(2)) resolution:resolution:mx(2)];
else
    grid.ygrid = mn(2):resolution:mx(2);
end

if mn(3)<0 && mx(3)>0
    grid.zgrid = [fliplr(0:-resolution:mn(3)) resolution:resolution:mx(3)];
else
    grid.zgrid = mn(3):resolution:mx(3);
end

if strcmp(side,'R')
   grid.xgrid = grid.xgrid(grid.xgrid>=0);
elseif strcmp(side,'L')
   grid.xgrid = grid.xgrid(grid.xgrid<=0);
end

grid.dim   = [length(grid.xgrid) length(grid.ygrid) length(grid.zgrid)];
[X, Y, Z]  = ndgrid(grid.xgrid, grid.ygrid, grid.zgrid);

pos   = [X(:) Y(:) Z(:)];

inside = ft_inside_headmodel(pos, struct('bnd', constraint));

if visualise_sources
    figure;
    title(['viewed in ' grid_space])
    ft_plot_mesh(constraint,'facealpha',0);hold on;
    p = pos(inside,:);
    plot3(p(:, 1),p(:, 2),p(:, 3), '.b', 'MarkerSize', 15);hold on;
end

grid.pos     = pos;
grid.inside  = find(inside);
grid.outside = find(~inside);
%grid.pos     = pos(inside, :);

%% output image setup
switch image_space
    case 'MNI template'
        sMRI   = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    case 'MNI-aligned'
        sMRI   = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    case 'native'
        sMRI   = D.inv{inv_index}.mesh.sMRI;
end
outvol       = spm_vol(sMRI);
outvol.dt(1) = spm_type('float32');
for i = 1:numel(image_name)
    outvol.fname = image_name{i};
    ov{i}        = spm_create_vol(outvol);
end


%% start the beamformer
if ROI
    % output variable
    out = [];
    parfor (n = 1:numel(inside),6)

        %spm('defaults','eeg');
        if ~inside(n)
            continue
        end
        disp(['gridpoint '  num2str(n) ' at co-ordinates ' num2str(pos(n,:))])


        %% prepare ROI data

        trans           = D.inv{1}.forward.fromMNI;
        head_roi_ctf    = spm_eeg_inv_transform_points(trans,pos(n,:));  % CTF co-ordinates
        fid             = D.fiducials;
        fid             = ft_convert_units(fid, 'm');

        %Mc = spm_eeg_inv_headcoordinates(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));

        M = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
        head_roi_elekta = spm_eeg_inv_transform_points(M,head_roi_ctf); % Centre of harmonic expansion

        S                = [];
        S.D              = D;
        S.r_sphere       = head_roi_elekta';
        S.Lin            = 8;
        S.roi_rad        = 0.03;
        S.corr_limit     = 0.98;
        S.input          = 'inner';
        S.roi_type       = 'sphere';
        S.prefix         = 'ROI_sphere_';
        S.interference_space = 'extern_';
        S.temporal_samples_threshold = 2.5e4;
        S.apply_temporal = 1;
        % new option to normalise weights
        S.normalise      = 1;
        [D_roi,~,weight]     = tsss_roi_ext(S);

        %% custom - optional epoching code

        D_roi                = custom_epochs(D_roi,stimfrq,peak_offset);


        %%
        matlabbatch = {};
        matlabbatch{1}.spm.tools.beamforming.data.dir = {bf_dir};
        matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(D_roi)};
        matlabbatch{1}.spm.tools.beamforming.data.val = 1;
        matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
        matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
        matlabbatch{1}.spm.tools.beamforming.data.overwrite = 1;
        matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
        matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
        matlabbatch{2}.spm.tools.beamforming.sources.plugin.mni_coords.pos = pos(n,:);
        matlabbatch{2}.spm.tools.beamforming.sources.normalise_lf = false;
        matlabbatch{2}.spm.tools.beamforming.sources.visualise = true;
        matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
        matlabbatch{3}.spm.tools.beamforming.features.woi = timewin;
        matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEG'};
        matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
        matlabbatch{3}.spm.tools.beamforming.features.cross_terms = 'megeeg';
        matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.foi = freqoi;
        matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.taper = 'hanning';
        matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 1;
        matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
        matlabbatch{3}.spm.tools.beamforming.features.visualise = 0;
        matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
        matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
        matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
        matlabbatch{5}.spm.tools.beamforming.output.plugin.montage.method = 'max';
        matlabbatch{5}.spm.tools.beamforming.output.plugin.montage.vois = cell(1, 0);
        matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));

        % matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'write';
        % matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'MEG';
        % matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.none = 0;
        % matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.prefix = 'B';

        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'online';
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'MEG';
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.none = 0;
        matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.prefix = 'B';
        spm_jobman('run',matlabbatch);

        % [pre,post] = spm_fileparts(spmfile);
        % D          = spm_eeg_load(fullfile(pre,['B' post]));
        % D          = copy(D,['C:\home\Data\4AshO\ROI_images\' 'B' num2str(n) post '.mat']);

        D_roi    = spm_eeg_load(fullfile(D_roi));
        mont     = montage(D_roi,'getmontage',montage(D_roi,'getnumber'));
        mont.tra = mont.tra./norm(mont.tra);
        D_roi    = montage(D_roi,'remove',montage(D_roi,'getnumber'));
        D_roi    = montage(D_roi,'add',mont);

        % compute power
        if strcmp(output.type,'spectra')
            data           = ftraw(D_roi,':',find(D_roi.time>output.time(1)*1e-3 & D_roi.time<output.time(2)*1e-3),':');
            cfg            = [];
            cfg.output     ='pow';
            cfg.keeptrials = 'no';
            cfg.keeptapers ='no';
            cfg.taper      = 'dpss';
            cfg.method     = 'mtmfft';
            cfg.foi        = freqoi(1):1:freqoi(2);
            cfg.tapsmofrq  = 5;
            cfg.pad        ='nextpow2';
            inp            = ft_freqanalysis(cfg, data);

            find = inp.freq>=frqout(1) & inp.freq<frqout(2);
            %oo = mean(inp.powspctrm(:,find),2);
            out = [out; mean(inp.powspctrm(:,find),2)];
        elseif strcmp(output.type,'rms')
            dat = D_roi(:,D_roi.time>output.time(1) & D_roi.time<outout.time(2),:);
            %oo = sqrt(mean(dat.^2,[2,3]));
            out = [out; sqrt(mean(dat.^2,[2,3]))];
        elseif strcmp(output.type,'peakamp')
            oo  = nan(1,numel(output.time))
            for ot = 1:numel(output.time)
                [~,i] = min(abs(D_roi.time*1e3-output.time(ot)));
                oo(ot) = mean(D_roi(:,i(1),:),3);
            end
            out = [out; oo];

        end

        % merge
        % if n<=1
        %     names  =[];
        % end
        % names = strvcat(names,fullfile(D));
        % if n>1
        %     S.D      = names;
        %     S.prefix = '';
        %     D        = spm_eeg_fuse(S);
        %     names    = fullfile(D);
        %     %delete(spm_eeg_load(S.D(1,:)));
        %     delete(spm_eeg_load(S.D(2,:)));
        % end

        % cleanup
        delete *BF.mat;
        %clear matlabbatch;
        delete(D_roi);

    end
else % do lcmv
    matlabbatch{1}.spm.tools.beamforming.data.dir = {bf_dir};
    matlabbatch{1}.spm.tools.beamforming.data.D = {spmfile};
    matlabbatch{1}.spm.tools.beamforming.data.val = 1;
    matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
    matlabbatch{1}.spm.tools.beamforming.data.space = 'MNI-aligned';
    matlabbatch{1}.spm.tools.beamforming.data.overwrite = 1;
    matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
    matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
    matlabbatch{2}.spm.tools.beamforming.sources.plugin.mni_coords.pos = grid.pos(inside,:);
    matlabbatch{2}.spm.tools.beamforming.sources.normalise_lf = false;
    matlabbatch{2}.spm.tools.beamforming.sources.visualise = true;
    matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
    matlabbatch{3}.spm.tools.beamforming.features.woi = timewin;
    matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEG'};
    matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
    matlabbatch{3}.spm.tools.beamforming.features.cross_terms = 'megeeg';
    matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.foi = freqoi;
    matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.taper = 'hanning';
    matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 1;
    matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
    matlabbatch{3}.spm.tools.beamforming.features.visualise = 0;
    matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
    matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
    matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{5}.spm.tools.beamforming.output.plugin.montage.method = 'max';
    matlabbatch{5}.spm.tools.beamforming.output.plugin.montage.vois = cell(1, 0);
    matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
    matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.mode = 'online';
    matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.modality = 'MEG';
    matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.addchannels.none = 0;
    matlabbatch{6}.spm.tools.beamforming.write.plugin.spmeeg.prefix = 'B';
    spm_jobman('run',matlabbatch);


    D    = spm_eeg_load(spmfile);
    mon_num = montage(D,'getnumber');
    mont = montage(D,'getmontage',mon_num);
    for i = 1:size(mont.tra,1)
        mont.tra(i,:) = mont.tra(i,:)./norm(mont.tra(i,:));
    end
    D    = montage(D,'remove',mon_num);
    D    = montage(D,'add',mont);
    D    = montage(D,'switch', mon_num);

    % compute power
    if strcmp(output.type,'spectra')
        data           = ftraw(D,':',find(D.time>output.time(1)*1e-3 & D.time<output.time(2)*1e-3),':');
        cfg            = [];
        cfg.output     ='pow';
        cfg.keeptrials = 'no';
        cfg.keeptapers ='no';
        cfg.taper      = 'dpss';
        cfg.method     = 'mtmfft';
        cfg.foi        = freqoi(1):1:freqoi(2);
        cfg.tapsmofrq  = 5;
        cfg.pad        ='nextpow2';
        inp            = ft_freqanalysis(cfg, data);

        find = inp.freq>=frqout(1) & inp.freq<frqout(2);
        out = mean(inp.powspctrm(:,find),2);

    elseif strcmp(output.type,'rms')

        dat = D(:,D.time>output.time(1) & D.time<outout.time(2),:);
        out = sqrt(mean(dat.^2,[2,3]));

    elseif strcmp(output.type,'peakamp')
        out = [];
        for ot = 1:numel(output.time)
            [~,i] = min(abs(D.time*1e3-output.time(ot)));
            out   = [out mean(D(:,i(1),:),3)];
        end

    end
    %D    = montage(D,'remove',mon_num);
    D    = montage(D,'switch',0);
    save(D);
end

%%

switch image_space
    % MNI template required no transform
    case 'MNI-aligned'
        grid = ft_transform_geometry(Data.transforms.toMNI_aligned/Data.transforms.toMNI, grid);
    case 'native'
        grid = ft_transform_geometry(BF.data.transforms.toNative/Data.transforms.toMNI, grid);
end

% separate image for each column
if normpow
    out = out./mean(out,1);
end

cfg                       = [];
cfg.parameter             = 'pow';
cfg.downsample            = 1;
cfg.showcallinfo          = 'no';
for o = 1:size(out,2)
grid.pow                  = nan(size(grid.pos, 1), 1);
grid.pow(grid.inside)     = out(:,o);
sourceint                 = ft_sourceinterpolate(cfg, grid, ft_read_mri(sMRI, 'dataformat', 'nifti_spm'));
Y                         = sourceint.pow;
spm_write_vol(ov{o}, Y);
end


%% custom epoching

function D = custom_epochs(D,stimfrq,peak_offset)

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

S           = [];
S.D         = D;
S.trl       = trl;
S.bc        = 0;
D           = spm_eeg_epochs(S);
D           = timeonset(D,-peak_offset*1e-3);
save(D);

if ~isfield(D,'inv')
    warning('running headmodel as this has not been specified');
    D = dbs_meg_headmodelling(D);
end
end
