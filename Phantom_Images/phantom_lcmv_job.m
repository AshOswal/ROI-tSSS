
function phantom_lcmv_job

MNI_roi    = [40 -20 50]; % this is the simulated dipole of interest

D          = spm_eeg_load('D:\home\Data\4AshO\Phantom Images\ROI_sim.mat');

%%
S            = [];
S.D          = D;
S.tsss       = 0;
S.t_window   = 5;
S.xspace     = 0;
D2           = tsss_spm_enm(S);

%%
trans           = D.inv{1}.forward.fromMNI;
head_roi_ctf    = spm_eeg_inv_transform_points(trans,MNI_roi);  % CTF co-ordinates
fid             = D.fiducials;
fid             = ft_convert_units(fid, 'm');

%Mc = spm_eeg_inv_headcoordinates(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));

M = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
head_roi_elekta = spm_eeg_inv_transform_points(M,head_roi_ctf); % Centre of harmonic expansion

S                = [];
S.D              = D2;
S.r_sphere       = head_roi_elekta';
S.Lin            = 8;
S.roi_rad        = 0.05;
S.corr_limit     = 1;
S.input          = 'inner';
S.roi_type       = 'sphere';
S.prefix         = 'ROI_tsss_sphere_';
S.interference_space = 'extern_';
S.temporal_samples_threshold = 2.5e4;
S.apply_temporal = 1;
% new option to normalise weights
S.normalise      = 1;
D3               = tsss_roi_ext(S);
D3 = forward_model_phantom(D3);
%%    

data = {D,D2,D3};

for d = 1:numel(data)
matlabbatch{1}.spm.tools.beamforming.data.dir = {'D:\home\Data\4AshO\Phantom Images'};
matlabbatch{1}.spm.tools.beamforming.data.D = {fullfile(data{d})};
matlabbatch{1}.spm.tools.beamforming.data.val = 1;
matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{1}.spm.tools.beamforming.data.space = 'Head';
matlabbatch{1}.spm.tools.beamforming.data.overwrite = 1;
matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.grid_phantom.resolution = 5;
matlabbatch{2}.spm.tools.beamforming.sources.normalise_lf = false;
matlabbatch{2}.spm.tools.beamforming.sources.visualise = 1;
matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{3}.spm.tools.beamforming.features.woi = [-Inf Inf];
matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEG'};
matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{3}.spm.tools.beamforming.features.cross_terms = 'megeeg';
matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.foi = [10 40];
matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.taper = 'none';
matlabbatch{3}.spm.tools.beamforming.features.regularisation.manual.lambda = 5;
matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{3}.spm.tools.beamforming.features.visualise = 1;
matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.orient = true;
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.lcmv.keeplf = false;
matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.whatconditions.all = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.sametrials = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.woi = [-Inf Inf];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.foi = [10 40];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.contrast = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.logpower = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.result = 'singleimage';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.scale = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.powermethod = 'trace';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_power.modality = 'MEG';
matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.write.plugin.nifti.normalise = 'no';
matlabbatch{6}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';
spm_jobman('run',matlabbatch);
end


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

