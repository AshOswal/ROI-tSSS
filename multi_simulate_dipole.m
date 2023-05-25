
clear all;
close all;

spmfile = 'D:\home\Data\DBS-MEG\Phantom\090715\SPMstim\spmeeg_uespmeeg_phantom090715_BrainampDBS_20150709_02.mat';
%spmfile = 'D:\home\Data\DBS-MEG\LN_C54\100316\SPMrest\LN_C54_on.mat';

D = spm_eeg_load(spmfile);

trial_length   = 40;      % (s)
ntrials        = 1;       % number of trials
dipfrq         = 20;      % dipole frequency
add_lfp        = 1;


sample_frq     = 200;     % new data sample frequency
time           = 0:1/sample_frq:trial_length;
dipamp         = 8e-2;    % fT
megsnr         = 2;

Data           = spm_eeg_inv_get_vol_sens(D,1,'Head','inv','MEG');

grad           = Data.MEG.sens;
vol            = Data.MEG.vol;
datareg        = D.inv{D.val}.datareg(1);

ROI_size        = 0.02; 
dipolepositions = [-0.03 0.02 0.03;0.03 0.02 0.03];

[x,y,z]         = meshgrid(-0.12:0.02:0.12,-0.12:0.02:0.12,-0.12:0.02:0.12);
noisepositions  = [x(:), y(:), z(:)];
%noisepositions  = [];

for n         = 1:size(dipolepositions,1)
oris          = [0 0 1]; % source orientation
oris          = oris/norm(oris);
oris          = spm_orth([dipolepositions(n,:)',oris']);
oris          = oris(:,2);
all_oris(n,:) = oris/norm(oris);
end

%%
%the skull mesh - transform the geometry
iskull = export(gifti(D.inv{D.val}.mesh.tess_iskull), 'ft');
trans = Data.transforms.toHead/Data.transforms.toNative;
iskull = ft_convert_units(ft_transform_geometry(trans, iskull),'m');

% make a plot
inside = ft_inside_headmodel(dipolepositions, struct('bnd', iskull));
dipolepositions = dipolepositions(find(inside),:);
inside = ft_inside_headmodel(noisepositions, struct('bnd', iskull));
noisepositions  = noisepositions(find(inside),:);

Ndips=size(dipolepositions,1);

plot_oris = Data.MEG.vol.r;
origin    = Data.MEG.vol.o;

figure;
ft_plot_headmodel(vol, 'edgecolor', [0 0 0], 'facealpha', 0);hold on;
for i = 1:size(dipolepositions,1)    
hm.o = dipolepositions(i,:);
hm.r = 0.02;
try
noisepositions = noisepositions(~ft_inside_vol(noisepositions,hm),:);
end
ft_plot_headmodel(hm,'edgecolor', [1 0 0], 'facealpha',0);hold on;
end

scale = 0.6;
quiver3(dipolepositions(:,1), dipolepositions(:,2), dipolepositions(:,3), all_oris(:,1), all_oris(:,2), all_oris(:,3),scale,'k','LineWidth',3);hold on;
plot3(origin(1),origin(2),origin(3),'.b','MarkerSize',40);hold on;
plot3(dipolepositions(:, 1),dipolepositions(:, 2),dipolepositions(:, 3), '.b', 'MarkerSize', 40);hold on;
try
plot3(noisepositions(:, 1),noisepositions(:, 2),noisepositions(:, 3), '.g', 'MarkerSize', 40);hold on;
end

scale = 0;
quiver3(repmat(origin(1),Ndips,1),repmat(origin(2),Ndips,1),repmat(origin(3),Ndips,1),dipolepositions(:,1),dipolepositions(:,2),dipolepositions(:,3),scale,'.b')
rotate3d off;
axis off
axis vis3d
axis equal
view([0 0])
%%
dipolepositions = [dipolepositions;noisepositions];
for n = 1:size(dipolepositions,1)
    cfg             = [];
    cfg.headmodel   = vol;
    cfg.grad        = grad;
    cfg.dip.unit    = 'nA*m';
    cfg.chanunit    = repmat('fT', 1, 274);
    cfg.dip.pos     = dipolepositions(n,:);
    cfg.ntrials     = ntrials;
    
    
    % simulate sinusoids
    %w               = (n-1) * 2*pi/size(dipolepositions,1);
    w               = [0 -pi/4 -pi/3 -pi/2];
    
    if n <= 2
        meg_signal(n,:) = dipamp.*sin(2*pi*dipfrq*time + w(n)); % variable phase offset        
        cfg.dip.mom     = all_oris(n,:);
    else
        rmsq            = dipamp/(sqrt(2)*15);       % formula for RMS sinusoid
        meg_signal(n,:) = rmsq*randn(1,numel(time)); % noise with amplitude 1/20 of simulated sinusoid
        rori            = rand(3,1);
        rori            = rori./norm(rori);
        cfg.dip.mom     = rori;

    end

    
    for i=1:cfg.ntrials
        cfg.dip.signal{i}=meg_signal(n,:);
    end
    
    
    cfg.triallength = time;% seconds
    cfg.fsample     = sample_frq;% Hz
    onesampletime   = 1/sample_frq;
    cfg.reducerank  = 1;
    cfg.absnoise    = 0;
    cfg.relnoise    = 0; % an arbitrary relative noise factor
    cfg.channel     = D.chanlabels(D.indchantype('MEG'));
    raws{n}         = ft_dipolesimulation(cfg);
    
    %%
    if ~isempty(megsnr)

    rms_raws_sq     = mean(mean(mean(cat(3,raws{n}.trial{:}).^2,3),2),1);
    noise_var       = rms_raws_sq/megsnr;
    
    % only add noise for signal of interest and not for noisy dipoles
        for i=1:cfg.ntrials
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
        raw.trial{i}   = [raw.trial{i}; cfg.dip.signal{:}];
        raw.label{275}   = 'LFP_CLEAN';
    end
end

% now write out a new data set
D2=spm_eeg_ft2spm(raw,'ROI_sim.mat');
D2=sensors(D2,'MEG',D.sensors('MEG'));
D2=fiducials(D2,D.fiducials);
D2.inv=D.inv;
D2=coor2D(D2,'MEG',coor2D(D));
D2 = units(D2, D2.indchantype('MEG'), 'fT');
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
cfg=[];
cfg.appenddim        = 'chan';
cfg.headmodel        = vol;
cfg.grad             = grad;
cfg.grid.resolution  = 0.01;
cfg.gridsearch       = 'yes';
cfg.channel          = 'all';
cfg.model            = 'regional';
cfg.reducerank       = 2;
cfg.numdipoles       = 2;%size(dipolepositions,1);
cfg.dipfit.noisecov  = cov(D2n(1:274,:)');
cfg.symmetry         = 'x';
dip_source           = ft_dipolefitting(cfg, d2);

idx       = find(dip_source.dip.pos(:,1)<0);

[u, s, v]    = svd(dip_source.dip.mom(1:3,:));
[u1, s1, v1] = svd(dip_source.dip.mom(4:6,:));

if idx==1
est_ori_d1   = u(:,1)'; % orientation
signal_d1    = s(1,1)*v(:,1)';
est_ori_d2   = u1(:,1)'; % orientation
signal_d2    = s1(1,1)*v1(:,1)';
else
est_ori_d1   = u1(:,1)'; % orientation
signal_d1    = s1(1,1)*v1(:,1)';
est_ori_d2   = u(:,1)'; % orientation
signal_d2    = s(1,1)*v(:,1)';    
end
% if angle between estimated vectors is >90 then flip 180 degrees
error_ori_d1 = atan2(norm(cross(est_ori_d1, all_oris(1,:))),dot(est_ori_d1, all_oris(1,:)));
error_ori_d2 = atan2(norm(cross(est_ori_d2, all_oris(2,:))),dot(est_ori_d2, all_oris(2,:)));
% if abs(error_ori_d1)>pi/2
%    keyboard;
%    est_ori_d1 = -est_ori_d1;
%    signal_d1  = -signal_d1;
% end
% if abs(error_ori_d2)>pi/2
%     keyboard;
%    est_ori_d2 = -est_ori_d2;
%    signal_d2  = -signal_d2;
% end

figure;
plot(meg_signal(1,1:100),'r');hold on;
plot(meg_signal(2,1:100),'b');hold on;
plot(signal_d1(1:100),'k');hold on;
plot(signal_d2(1:100),'g');hold on
legend({'true dipole 1','true dipole 2','estimated dipole 1','estimated dipole 2'});

% error_ori = min(abs(error_ori), abs(error_ori-pi/2));
% error_pos = norm(dip_source.dip.pos - dipolepositions);
%%
figure;
ft_plot_dipole(dipolepositions(1,:), all_oris(1,:), 'color', 'r','unit','m','alpha',1 );hold on
ft_plot_dipole(dipolepositions(2,:), all_oris(2,:), 'color', 'r','unit','m','alpha',1  );hold on
ft_plot_dipole(dip_source.dip.pos(1,:), est_ori_d1(1,:), 'color', 'g','unit','m','alpha',1 );hold on
ft_plot_dipole(dip_source.dip.pos(2,:), est_ori_d2(1,:), 'color', 'g','unit','m','alpha',1 );hold on
rotate3d off;
axis off
axis vis3d
axis equal
view([0 0])

% if nl ==0
%     save('est_pos','est_pos');
%     save('est_ori','est_ori');
% end
