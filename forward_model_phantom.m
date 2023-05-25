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