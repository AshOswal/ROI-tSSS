%-----------------------------------------------------------------------
% Job saved on 19-May-2023 21:06:31 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.beamforming.data.dir = {'D:\home\Data\4AshO'};
matlabbatch{1}.spm.tools.beamforming.data.D = {'D:\home\Data\DBS-MEG\Phantom\090715\SPMstim\efsss_spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_11.mat'};
 
% 'efROI_tsss_sphere_einner_rsss_spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_11.mat']};

matlabbatch{1}.spm.tools.beamforming.data.val = 1;
matlabbatch{1}.spm.tools.beamforming.data.gradsource = 'inv';
matlabbatch{1}.spm.tools.beamforming.data.space = 'Head';
matlabbatch{1}.spm.tools.beamforming.data.overwrite = 0;
matlabbatch{2}.spm.tools.beamforming.sources.BF(1) = cfg_dep('Prepare data: BF.mat file', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{2}.spm.tools.beamforming.sources.reduce_rank = [2 3];
matlabbatch{2}.spm.tools.beamforming.sources.keep3d = 1;
matlabbatch{2}.spm.tools.beamforming.sources.plugin.grid_phantom.resolution = 2;
matlabbatch{2}.spm.tools.beamforming.sources.normalise_lf = false;
matlabbatch{2}.spm.tools.beamforming.sources.visualise = 1;
matlabbatch{3}.spm.tools.beamforming.features.BF(1) = cfg_dep('Define sources: BF.mat file', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{3}.spm.tools.beamforming.features.whatconditions.all = 1;
matlabbatch{3}.spm.tools.beamforming.features.woi = [-Inf Inf];
matlabbatch{3}.spm.tools.beamforming.features.modality = {'MEG'};
matlabbatch{3}.spm.tools.beamforming.features.fuse = 'no';
matlabbatch{3}.spm.tools.beamforming.features.cross_terms = 'megeeg';
matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.foi = [20 30];
matlabbatch{3}.spm.tools.beamforming.features.plugin.cov.taper = 'hanning';
matlabbatch{3}.spm.tools.beamforming.features.regularisation.mantrunc.pcadim = 75;
matlabbatch{3}.spm.tools.beamforming.features.bootstrap = false;
matlabbatch{3}.spm.tools.beamforming.features.visualise = 1;
matlabbatch{4}.spm.tools.beamforming.inverse.BF(1) = cfg_dep('Covariance features: BF.mat file', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{4}.spm.tools.beamforming.inverse.plugin.dics.fixedori = 'yes';
matlabbatch{5}.spm.tools.beamforming.output.BF(1) = cfg_dep('Inverse solution: BF.mat file', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.reference.power = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.powmethod = 'lambda1';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.whatconditions.all = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.sametrials = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.woi = [-Inf Inf];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.contrast = 1;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.logpower = false;
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.foi = [20 30];
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.taper = 'hanning';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.result = 'singleimage';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.scale = 'yes';
matlabbatch{5}.spm.tools.beamforming.output.plugin.image_dics.modality = 'MEG';
matlabbatch{6}.spm.tools.beamforming.write.BF(1) = cfg_dep('Output: BF.mat file', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','BF'));
matlabbatch{6}.spm.tools.beamforming.write.plugin.nifti.normalise = 'separate';
matlabbatch{6}.spm.tools.beamforming.write.plugin.nifti.space = 'mni';
spm_jobman('run',matlabbatch)
