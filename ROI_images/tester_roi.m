addpath('C:\home\Data\4AshO')
D = spm_eeg_load()
MNI_roi = [40 -20 40];

trans           = D.inv{1}.forward.fromMNI;
head_roi_ctf    = spm_eeg_inv_transform_points(trans,MNI_roi);  % CTF co-ordinates
fid             = D.fiducials;
fid             = ft_convert_units(fid, 'm');

%Mc = spm_eeg_inv_headcoordinates(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));

M = headcoordinates_mne(fid.fid.pnt(1, :), fid.fid.pnt(2, :), fid.fid.pnt(3, :));
head_roi_elekta = spm_eeg_inv_transform_points(M,head_roi_ctf); % Centre of harmonic expansion

S            = [];
S.D          = D;
S.tsss       = 1;
S.t_window   = 5;
S.xspace     = 0;
Dt           = tsss_spm_enm(S);
Dt           = dbs_meg_headmodelling(Dt);
Dt.save;

S                = [];
S.D              = Dt;
S.r_sphere       = head_roi_elekta';
S.Lin            = 8;
S.roi_rad        = 0.05;
S.corr_limit     = 1;
S.input          = 'inner';
S.roi_type       = 'sphere';
S.prefix         = 't_';
S.interference_space = 'outer';
S.temporal_samples_threshold = 2.5e4;
S.apply_temporal = 1;
Do               = tsss_roi_ext(S);
Do               = dbs_meg_headmodelling(Do);
Do.save;