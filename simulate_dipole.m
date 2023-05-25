
clear all;
close all;
%spm('defaults','eeg');
spmfile = 'D:\home\Data\DBS-MEG\Phantom\090715\SPMstim\spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_02.mat';
%spmfile = 'D:\home\Data\DBS-MEG\LN_C54\100316\SPMrest\LN_C54_on.mat';
%spmfile = 'D:\home\Data\DBS-MEG\SW\190207\SPMrest\SW_off.mat';
D = spm_eeg_load(spmfile);

trial_length   = 40;         % (s)
ntrials        = 1;           % number of trials
dipoleposition = [40 -20 40;];% 20 -20 40;-40 -20 40;50 20 30;-50 20 30; -20 -20 40]; % dipole location in MNI space
dipfrq         = [20];%;15;19;25;28;30];        % dipole frequncy
add_lfp        = 1;


sample_frq     = 200;     % new data sample frequency
time           = 0:1/sample_frq:trial_length;
dipamp         = 8e-2;    % fT
megsnr         = 0.5;

Data           = spm_eeg_inv_get_vol_sens(D,1,'MNI-aligned','inv','MEG');

grad           = Data.MEG.sens;
vol            = Data.MEG.vol;
datareg        = D.inv{D.val}.datareg(1);

M = Data.transforms.toMNI;
dipolepositions =  spm_eeg_inv_transform_points(inv(M),dipoleposition);

oris          = [0 0 1]; % source orientation
oris          = oris/norm(oris);
oris          = spm_orth([dipolepositions(1,:)',oris']);
oris          = oris(:,2);
oris          = oris/norm(oris);

%%
%the skull mesh - transform the geometry
iskull = export(gifti(D.inv{D.val}.mesh.tess_iskull), 'ft');
trans = Data.transforms.toMNI_aligned/Data.transforms.toNative; % affine transform from Native to MNI_aligned
iskull = ft_convert_units(ft_transform_geometry(trans, iskull),'m');

% make a plot
inside = ft_inside_headmodel(dipolepositions, struct('bnd', iskull));
dipolepositions = dipolepositions(find(inside),:);
Ndips=size(dipolepositions,1);


figure;
ft_plot_headmodel(vol, 'edgecolor', [0 0 0], 'facealpha', 0);hold on;
plot3(dipolepositions(:, 1),dipolepositions(:, 2),dipolepositions(:, 3), '.r', 'MarkerSize', 30);hold on;
rotate3d off;
axis off
axis vis3d
axis equal
view([270 0])



%%

for n = 1:size(dipolepositions,1)
    cfg             = [];
    cfg.headmodel   = vol;
    cfg.grad        = grad;
    cfg.dipoleunit  = 'nA*m';
    cfg.chanunit    = 'fT';%repmat('fT', 1, 274);
    cfg.sourcemodel.pos     = dipolepositions(n,:);
    cfg.numtrl      = ntrials;
    
    
    % simulate sinusoids
    meg_signal      = dipamp*sin(time*dipfrq(n)*2*pi);
  
    %cfg.sourcemodel.mom = oris;
    cfg.sourcemodel.mom = [0 1 0];
    
    for i=1:cfg.numtrl
        cfg.sourcemodel.signal{i}=meg_signal;
    end
    
    
    cfg.trllen      = time;% seconds
    cfg.fsample     = sample_frq;% Hz
    onesampletime   = 1/sample_frq;
    cfg.reducerank  = 2;
    cfg.absnoise    = 0;
    cfg.relnoise    = 0; % an arbitrary relative noise factor
    cfg.channel     = D.chanlabels(D.indchantype('MEG'));
    raws{n}         = ft_dipolesimulation(cfg);
    
    if ~isempty(megsnr)

    rms_raws_sq     = mean(mean(mean(cat(3,raws{n}.trial{:}).^2,3),2),1);
    noise_var       = rms_raws_sq/megsnr;
    
    % only add noise for signal of interest and not for noisy dipoles
        for i=1:cfg.numtrl
            if isequal(n,1)
                noise  = sqrt(noise_var)*randn(size(raws{n}.trial{i}));
            else
                noise  = zeros(size(raws{n}.trial{i}));
            end
            raws{n}.trial{i}    = raws{n}.trial{i} + noise;
            raws{n}.noise{i}    = noise;
            raws{n}.noisecov{i} = cov(noise');
        end
    end
    
    rn = cfg.relnoise;
    an = cfg.absnoise;
    if isequal(an,0) && isequal(rn,0)
        nl = 0;
    else
        nl = [];
    end
end
% merge into one
raw = raws{1};
for u = 1:numel(raw.trial)
    raw.trial{u} = zeros(size(raw.trial{u}));
    if ~isempty(megsnr)
    raw.noise{u} = zeros(size(raw.noise{u}));
    end
end

for n = 1:numel(raws)
    for t = 1:numel(raw.trial)
        raw.trial{t} =  raw.trial{t} + raws{n}.trial{t};
        if ~isempty(megsnr)
        raw.noise{t} =  raw.noise{t} + raws{n}.noise{t};
        end
    end
end

if  add_lfp
    for i = 1:numel(raw.trial)
        raw.trial{i}   = [raw.trial{i}; cfg.sourcemodel.signal{:}];
        raw.label{end+1}   = 'LFP_CLEAN';
    end
end

% now write out a new data set
D2=spm_eeg_ft2spm(raw,'ROI_sim.mat');
D2=sensors(D2,'MEG',D.sensors('MEG'));
D2=fiducials(D2,D.fiducials);
D2.inv=D.inv;
D2=coor2D(D2,'MEG',coor2D(D));
D2 = units(D2, D2.indchantype('MEG'), 'fT');
if 1
    D2.initials = 'SW';
end
D2.save;

% separately write out a noise dataset 
if ~isempty(megsnr)
raw.trial      = raw.noise;
raw.label(275) = [];
D2n=spm_eeg_ft2spm(raw,'ROI_noise.mat');
D2n=sensors(D2n,'MEG',D.sensors('MEG'));
D2n=fiducials(D2n,D.fiducials);
D2n.inv=D.inv;
D2n=coor2D(D2n,'MEG',coor2D(D));
D2n = units(D2n, D2n.indchantype('MEG'), 'fT');
D2n.save;
end


%% epoching and dipole fit

d2   = ftraw(D2,setdiff(D2.indchantype('MEG'),D2.badchannels),':',':');
d2_s = ftraw(D2,setdiff(D2.indchantype('MEG'),D2.badchannels),3:4,':');

cfg            = [];
cfg.channel    = 'all';
cfg.trials     = 'all';
cfg.keeptrials = 'no';
d2_s             = ft_timelockanalysis(cfg,d2_s);

%%
cfg            = [];
cfg.parameter  = 'avg';
cfg.colorbar   = 'east';
cfg.style      = 'straight';
cfg.marker     = 'off';
%cfg.zlim       = [-5 5];
cfg.comment    = 'no';
figure;
ft_topoplotER(cfg,d2_s);

%% compare simulated and reconstructed timecourse

%data = spm_eeg_inv_get_vol_sens(D,1,'MNI-aligned','inv','MEG');

cfg=[];
cfg.appenddim  = 'chan';
cfg.headmodel  = vol;
cfg.grad = grad;
cfg.grid.resolution = 0.01;
cfg.gridsearch = 'yes';
cfg.channel    = 'all';
cfg.model      = 'regional';
cfg.reducerank = 2;
% various inverse options all work well
%cfg.weight     = inv(cov(D2n(1:274,:)'));

dip_source     = ft_dipolefitting(cfg, d2);

[u, s, v] = svd(dip_source.dip.mom);
est_ori = u(:,1)'; % orientation
signal = s(1,1)*v(:,1)';
%%

figure;plot(time,signal);hold on;plot(time,meg_signal,'r');
legend('reconstructed','simulated')

error_ori = atan2(norm(cross(est_ori, oris)),dot(est_ori, oris));

error_ori = min(abs(error_ori), abs(error_ori-pi/2));

error_pos = norm(dip_source.dip.pos - dipolepositions);
%%
figure;
ft_plot_dipole(dipolepositions, oris, 'color', 'r' );hold on
ft_plot_dipole(dip_source.dip.pos, est_ori, 'color', 'g' );hold on
rotate3d off;
axis off
axis vis3d
axis equal
view([270 0])
est_pos = dip_source.dip.pos;

if nl ==0
    save('est_pos','est_pos');
    save('est_ori','est_ori');
end
