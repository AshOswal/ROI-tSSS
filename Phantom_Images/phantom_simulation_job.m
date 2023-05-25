%-----------------------------------------------------------------------
% Job saved on 20-May-2023 23:13:37 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.beamforming.data.dir = {'D:\home\Data\4AshO\Phantom Images'};
matlabbatch{1}.spm.tools.beamforming.data.D = {'D:\home\Data\4AshO\Phantom Images\ROI_sim.mat'};
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
