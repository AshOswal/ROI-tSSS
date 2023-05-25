function [Di,Do,spatial_weight] = tsss_roi_ext(S)


if ~isfield(S, 'Dref'),                        S.Dref              = '';    end
if ~isfield(S, 'magscale'),                    S.magscale          = 100;   end
if ~isfield(S, 'Lin'),                         S.Lin               = 8;     end
if ~isfield(S, 'Lout'),                        S.Lout              = 3;     end
if ~isfield(S, 'prefix'),                      S.prefix            = 'ROI_'; end
if ~isfield(S, 'roi_rad'),                     S.roi_rad           = 0; end
if ~isfield(S, 'apply_temporal'),              S.apply_temporal    = 0; end 
if ~isfield(S, 'corr_limit'),                  S.corr_limit        = 0.98;  end
if ~isfield(S, 'keep_files'),                  S.keep_files        = 1;  end
if ~isfield(S, 'apply_temporal'),              S.apply_temporal    = 1;  end
if ~isfield(S, 'temporal_samples_threshold'),  S.temporal_samples_threshold    = 1e4;  end
if ~isfield(S, 'interference_space'),          S.interference_space   = 'external';  end
if ~isfield(S, 'normalise'),                   S.normalise         = 0;  end

%%
spm('FnBanner', mfilename);

%%
if ~isfield(S, 'roi_type'),         S.roi_type    = 'cube'; end

centre_ROI   = S.r_sphere;
D            = spm_eeg_load(S.D);
mag2SI       = 1e-15;
grad2SI      = 1e-12;
Lin          = S.Lin;
Lout         = S.Lout;
magscale     = S.magscale;
roirad       = S.roi_rad;
keep         = S.keep_files;
% temporal correlation
corr_limit   = S.corr_limit;

isneuromag = strncmpi(ft_senstype(D.chanlabels), 'neuromag', 7);

%%
meg_ch   = D.indchantype({'MEG', 'MEGPLANAR'});
other_ch = setdiff(1:D.nchannels, meg_ch);

if isneuromag
    mags    = strmatch('MEGMAG', D.chantype(meg_ch));
    grads   = strmatch('MEGPLANAR', D.chantype(meg_ch));
else
    mags     = 1:length(meg_ch);
    grads    = [];
end

goodchs = find(~D.badchannels(meg_ch));
coils = zeros(1, length(meg_ch));
coils(mags) = 1;

if isfield(D, 'allsens')
    allsens = D.allsens;
    allfid  = D.allfid;
else
    allsens = D.sensors('MEG');
    allfid  = D.fiducials;
end

nlocations = numel(allsens);

for i = 1:nlocations
        
    disp('Calculating the ROI SSS basis...');
    
    if isneuromag
        [R,EX,EY,EZ] = origheader_getpos(origheader{i}, 'head');
    else 
        [R,EX,EY,EZ] = ft_getpos(allsens(i), allfid(i), D.chanlabels(meg_ch));%
    end
    
    [centre,Radius] = minboundsphere(R'); % R : big_radius and 'center' is the original center
    R_prime = Radius + norm(centre_ROI - centre');

    if isneuromag
        Sin = Sin_vsh_fast(centre_ROI,R,EX,EY,EZ,coils,Lin);
    else
        Sin = Sin_vsh_ctf(centre_ROI,R,EX,EY,EZ,Lin);
    end
      
    %Sin(mags,:) = magscale*Sin(mags,:);
    
    for j = 1:size(Sin,2)
        SNin(:,j) = Sin(:,j)/norm(Sin(:,j));
    end    
    
    %% add outer in case it is needed
    if Lout>0
        if isneuromag
            Sout = Sout_vsh_fast(centre_ROI,R,EX,EY,EZ,coils,Lout);
        else
            Sout = Sout_vsh_ctf(centre_ROI,R,EX,EY,EZ,Lout);
        end
        %Sout(mags,:) = magscale*Sout(mags,:);
    else
        Sout  = [];
        SNout = [];
    end
   
    for j = 1:size(Sout,2)
        SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
    end
    
    SN   = SNin;
    pSN  = pinv(SN(goodchs,:));
    SNo  = SNout;
    pSNo = pinv(SNo(goodchs,:)); 
    SSS(i) = struct('Sin', Sin, 'SNin', SNin, 'SN', SN, 'pSN', pSN,'Sout',Sout,'SNout',SNout,'pSNout',pSNo);
end

if isfield(D, 'fileind')
    fileind = D.fileind;
else
    fileind = ones(1, D.ntrials);
end    
        
if strcmp(S.roi_type,'sphere')
    
    if S.normalise
        [deep_weight,sup_weight] = obtain_weights(roirad,R_prime,Lin);
    else
        [deep_weight,sup_weight] = obtain_weights_unnormalised(roirad,R_prime,Lin);
    end
    
    deep_weight = diag(deep_weight);
    sup_weight  = diag(sup_weight);
   
    % we break this into an inner, outer and external subspace
    if strcmp(S.interference_space,'external')
       sup_weight = diag(1-diag(deep_weight)); 
    end
    
elseif strcmp(S.roi_type,'cube')
    %deep_weight = GsGd(S.Lin,S.Lin,[roirad*2 R_prime],'inner');
    if S.normalise
        [G,sup_scale,deep_scale]  = GsGd(S.Lin,S.Lin,[roirad*2 R_prime],1);
        deep_weight = (G./G(1,1)) * deep_scale;
        sup_weight = inv(G);
        % really should use the first equation, the second is a numerical
        % hack to prevent weights above 1
        % sup_weight = (sup_weight./sup_weight(end,end)) * sup_scale;
        sup_weight = (sup_weight./max(sup_weight(:))) * sup_scale;
    else
        [G,~,deep_scale]  = GsGd(S.Lin,S.Lin,[roirad*2 R_prime],0);
        deep_weight = G*deep_scale;
        sup_weight  = inv(deep_weight);
    end
    
    if strcmp(S.interference_space,'external')
        sup_weight = diag(1-diag(deep_weight));
    end
end

%%
t = unique(fileind);
for i = 1:numel(unique(fileind))
    % use normalised version or not
    SNin = SSS(i).Sin;
    %SNin = SSS(i).Sin;
    SNout = SSS(i).Sout;
    
    D1  = badtrials(D,find(fileind~=t(i)),1);
    S1   = [];
    S1.D = D1; 
    D1 = spm_eeg_remove_bad_trials(S1);
    
    %% A few comments
    
    % note that the tsss_spm_enm function uses magnetometer and gradiometer
    % scaling prior to computing spherical harmonic basis sets. In that
    % function, the scaling is first applied to the data before processing
    % and then unapplied. This may not be necessary however.
    
    % In the below code note that SNin = diag(x)*Sin, where x is the
    % normaliser. Note also that pinv(SNin) = 1/diag(x)*pinv(Sin). The
    % normaliser and its inverse finally cancel out
    %%
    
    disp('computing spatial projector...')
    S1    = [];
    S1.D  = D1;
    S1.montage.tra = real(SNin*deep_weight*pinv(SNin));%/(mag2SI*magscale);
    S1.montage.labelorg = D.chanlabels(D.indchantype('MEG'));
    S1.montage.labelnew = S1.montage.labelorg;
    S1.keepsensors = 1;
    S1.mode = 'write';
    S1.prefix = 'inner_';
    Din{i} = spm_eeg_montage(S1);
    
    spatial_weight{i} = S1.montage.tra;
    
    if S.apply_temporal
        spm_progress_bar('Set', 0.5) 

        S1    = [];
        S1.D  = D1;
        S1.montage.tra = real(SNin*sup_weight*pinv(SNin));%/(mag2SI*magscale);
        S1.montage.labelorg = D.chanlabels(D.indchantype('MEG'));
        S1.montage.labelnew = S1.montage.labelorg;
        S1.keepsensors = 1;
        S1.mode = 'write';
        S1.prefix = 'outer_';
        Dout{i} = spm_eeg_montage(S1);

    end
    
    delete(D1);
end

if i>1 % combine
    S1 = [];
    S1.D = fname(Din{1});
    for f = 2:numel(Din)
        S1.D = strvcat(S1.D, fname(Din{f}));
    end
    S1.recode = 'same';
    Din = spm_eeg_merge(S1);
    
    if S.apply_temporal
        S1 = [];
        S1.D = fname(Dout{1});
        for f = 2:numel(Dout)
            S1.D = strvcat(S1.D, fname(Dout{f}));
        end
        S1.recode = 'same';
        Dout = spm_eeg_merge(S1);
    end
    
else
    Din    = Din{1};
    if S.apply_temporal
        Dout  = Dout{1};
    end
end

% an equivalent to tSSS for regions of interest
if S.apply_temporal
    
    
    if Din.ntrials ==1 && Din.nsamples> S.temporal_samples_threshold
        % crop the datasets if the temporal projector is likely to be large
        S1 = [];
        S1.bc = 0;
        S1.D = Din;
        S1.trialength = S.temporal_samples_threshold/Din.fsample*1e3;
        Din = spm_eeg_epochs(S1); % innner
        if ~keep, delete(S1.D); end

        
        S1.D = Dout;
        Dout = spm_eeg_epochs(S1); % outer
        if ~keep, delete(S1.D); end

        S1.D = D;
        De   = spm_eeg_epochs(S1);
 
    else
        De = D; 
    end
    

    Di = clone(Din, spm_file(fullfile(Din), 'prefix', S.prefix), [Din.nchannels Din.nsamples Din.ntrials], 0);
    Do = clone(Dout, spm_file(fullfile(Dout), 'prefix', S.prefix), [Dout.nchannels Din.nsamples Din.ntrials], 0);
    
    spm_progress_bar('Init', Din.ntrials,'Applying temporal projector');
    
    for trial = 1:Din.ntrials
        
        idata = Din(meg_ch,:,trial);
        odata = Dout(meg_ch,:,trial);
        % residual
        rdata = De(meg_ch,:,trial) - (idata+odata);
        
        Bin  = idata./norm(idata);
        Bout = odata./norm(odata);
        
        QA = orth(Bin');     % right singular value
        QB = orth(Bout');   % right singular value
        %[QA,ignore] = qr(Ein,0);
        %[QB,ignore] = qr(Eres,0);
        [U,SS,V] = svd(QA'*QB);
        %Up = QA*U;
        Vp = QB*V;
        s = diag(SS);
        inter_indices = find(s>=corr_limit);
        length(inter_indices);
        Eproj = Vp(:,inter_indices);
        
        P = Eproj*Eproj';
        Xp = (eye(size(P,1))-P);
        spm_progress_bar('Set', trial);
        
        %{
        if 0
            [spectrum1, ntaper, freqoi] = ft_specest_mtmfft(idata,Din.time,'taper','dpss','tapsmofrq',2.5,'freqoi',[0:2.5:90]);
            [spectrum2, ntaper, freqoi] = ft_specest_mtmfft(rdata,Din.time,'taper','dpss','tapsmofrq',2.5,'freqoi',[0:2.5:90]);
            [spectrum3, ntaper, freqoi] = ft_specest_mtmfft(Eproj',Din.time,'taper','dpss','tapsmofrq',2.5,'freqoi',[0:2.5:90]);
            [spectrum4, ntaper, freqoi] = ft_specest_mtmfft(Din(meg_ch,:,trial)*Xp,Din.time,'taper','dpss','tapsmofrq',2.5,'freqoi',[0:2.5:90]);
           
            pow = squeeze(mean(abs(spectrum1),1));
            pow2 = squeeze(mean(abs(spectrum2),1));
            pow3 = squeeze(mean(abs(spectrum3),1));
            pow4 = squeeze(mean(abs(spectrum4),1));
            
            figure;
            subplot(1,4,1);imagesc(freqoi,1:274,log(pow));
            subplot(1,4,2);imagesc(freqoi,1:274,log(pow2));
            subplot(1,4,3);imagesc(freqoi,1:length(inter_indices),log(pow3));
            subplot(1,4,4);imagesc(freqoi,1:274,log(pow4));
        end
        %}
        
        % inner temporal dataset
        Di(meg_ch,:,trial)   = Din(meg_ch,:,trial)*Xp;
        Di(other_ch,:,trial) = Din(other_ch,:,trial);
        % outer temporal dataset
        Do(meg_ch,:,trial)   = Dout(meg_ch,:,trial)*Xp;
        Do(other_ch,:,trial) = Dout(other_ch,:,trial);

    end
    
    if ~keep, delete(Din), delete(Dout); end
    
    % if D was not epoched then convert back to continuous dataset - which
    % may be cropped slightly depending on threshold trial length
    
    if D.ntrials ==1 && Di.ntrials>1
        warning('input and output data time dimension may differ due to epoching for temporal processing')
        
        data_in = reshape(Di(:,:,:),size(Di(:,:,:),1),[]);
        data_out= reshape(Do(:,:,:),size(Do(:,:,:),1),[]);
        
        Di = clone(Di, spm_file(fullfile(Di), 'prefix', []), [Di.nchannels Di.nsamples* Di.ntrials 1], 0);
        Di(:,:)=data_in;
        
        Do = clone(Do, spm_file(fullfile(Do), 'prefix', []), [Do.nchannels Do.nsamples* Do.ntrials 1], 0);
        Do(:,:)=data_out;
        
        Do = type(Do, 'continuous');
        Di = type(Di, 'continuous');
    end
    
    spm_progress_bar('Clear');
else
    % copy new dataset
    Din = copy(Din,[S.prefix Din.fname]);
    Di = Din;
    Do = [];
    return;
end
%-Save the M/EEG dataset
%--------------------------------------------------------------------------
Do = Do.history(mfilename, S);
Di = Di.history(mfilename, S);
save(Do); save(Di);



function [deep_weight,sup_weight] = obtain_weights(separating_radius,sensor_radius,Lin)

i =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );
Ls = [];
for L=1:Lin
    for M=-L:L

        deep_weight(i) = (sensor_radius^5 - separating_radius^5) .*  ( (separating_radius^(2*L-2)) ./ ...
            ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  );

        % normalise and ensure that tends to unity as r > R
        Ls = [Ls,L];
        i = i + 1;
    end
end
%sup_weight  = (1./deep_weight)./(1/deep_weight(end)) .* ((2*Lin)+3) ./ (2.*Ls+3);
deep_weight = deep_weight .* ((2*Ls)+3)/5;

sup_weight = (separating_radius.^(2*Lin - 2*Ls) .* (sensor_radius.^(2*Ls +3)) - separating_radius^(2*Lin+3)) ./...
    (sensor_radius^(2*Lin +3) - separating_radius^(2*Lin +3));

sup_weight = sup_weight .* ((2*Lin)+3) ./ (2.*Ls+3);

%% OLD
%{
function [deep_weight,sup_weight] = obtain_weights(separating_radius,sensor_radius,Lin)

i =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );
Ls = [];
for L=1:Lin
    for M=[-L:L]
        %deep_weight(i) = (sensor_radius^3 - separating_radius^3) .*  ( (separating_radius^(2*L)) ./ ...
        %    ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  ) .* (2*L+3) ./ 3 ;
        
        deep_weight(i) = (sensor_radius^3 - separating_radius^3) .*  ( (separating_radius^(2*L)) ./ ...
             ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  );
        
        % normalise and ensure that tends to unity as r > R 
        Ls = [Ls,L];
        %deep_weight(i) = deep_weight(i).* (2*L+3)./5;
        i = i + 1;
    end
end
sup_weight  = (1./deep_weight)./(1/deep_weight(end)) .* ((2*Lin)+3) ./ (2.*Ls+3);
deep_weight = deep_weight./deep_weight(1) .* ((2*Ls)+3)/5;

%sup_weight = 1./deep_weight;
%sup_weight = 1-deep_weight;
%}


function [deep_weight,sup_weight] = obtain_weights_unnormalised(separating_radius,sensor_radius,Lin)

i =1;
deep_weight = zeros(1,(Lin+1)^2 - 1 );

for L=1:Lin
    for M=[-L:L]
        deep_weight(i) = (sensor_radius^3 - separating_radius^3) .*  ( (separating_radius^(2*L)) ./ ...
            ((sensor_radius^(2*L+3)) - (separating_radius^(2*L+3)))  ) .* (2*L+3) ./ 3 ;
        i = i + 1;
    end
end
sup_weight = 1./deep_weight;



