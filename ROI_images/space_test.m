% dip test

% this code tests that the elekta co-ordinate system for SSS basis sets
% works...by imaging simulated sources within it

clear all;
close all;
spm('defaults','eeg');addpath 'C:\home\Data\4AshO'
%spmfile = 'C:\home\Data\DBS-MEG\Phantom\090715\SPMstim\spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_02.mat';
% spmfile = 'C:\home\Data\DBS-MEG\LN_C54\100316\SPMrest\LN_C54_on.mat';
spmfile = 'D:\home\Data\DBS-MEG\SW\190207\SPMrest\SW_off.mat';
D = spm_eeg_load(spmfile);

dipoleposition = [40 -20 40;0 -20 40;-40 -20 40;20 -20 40]; % dipole location in MNI space
Data           = spm_eeg_inv_get_vol_sens(D,1,'Head','inv','MEG');

grad           = Data.MEG.sens;
vol            = Data.MEG.vol;
datareg        = D.inv{D.val}.datareg(1);

M = Data.transforms.toMNI;
dipolepositions =  spm_eeg_inv_transform_points(inv(M),dipoleposition);

fid             = D.fiducials;
fid             = ft_convert_units(fid, 'm');
M               = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
dipolepositions = spm_eeg_inv_transform_points(M,dipolepositions); 

%%
%the skull mesh - transform the geometry
iskull = export(gifti(D.inv{D.val}.mesh.tess_iskull), 'ft');
trans = M/(Data.transforms.toHead/Data.transforms.toNative); % affine transform from Native to MNI_aligned
iskull = ft_convert_units(ft_transform_geometry(trans, iskull),'m');

% make a plot
inside = ft_inside_headmodel(dipolepositions, struct('bnd', iskull));
dipolepositions = dipolepositions(find(inside),:);
Ndips=size(dipolepositions,1);

% apply transform to volume conductor model
vol  = ft_transform_geometry(M, vol);

figure;
ft_plot_headmodel(vol, 'edgecolor', [0 0 0], 'facealpha', 0);hold on;
plot3(dipolepositions(:, 1),dipolepositions(:, 2),dipolepositions(:, 3), '.r', 'MarkerSize', 30);hold on;
rotate3d off;
axis off
axis vis3d
axis equal
view([270 0])

