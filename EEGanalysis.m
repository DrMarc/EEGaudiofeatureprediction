function EEGanalysis(eegfile,step)

% determine step
if nargin < 2
    step = 1;
    while exist(sprintf('.%s%s_step%i.set',filesep,eegfile(1:end-4),step),'file');
        step = step + 1;
    end
end
 
switch step
    case 1
        global EEG
        disp('Step 1: Manual pruning of continuous data');
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        EEG = pop_biosig(eegfile,'channels',[1:3,8:22]);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG.chanlocs = readlocs('chanlocs.ced'); 

        % rereference to linked mastoids
        EEG.data(16,:) = EEG.data(16,:)/2;
        EEG = pop_reref(EEG,16);
        EEG = pop_select(EEG,'nochannel',11); % rigt DRL

        % filter
        EEG = pop_eegfiltnew(EEG,[],1,8250,true,[],0);
        EEG = pop_eegfiltnew(EEG,[],40,330,0,[],0);
        EEG.setname='Step1result';

        % now reject by eye
        pop_eegplot(EEG,1,0,1);
        assignin('base','EEG',EEG); % move EEG variable to base workspace because pop_eegplot is not blocking
        uiwait;
        EEG = evalin('base','EEG'); % recover modified EEG variable from base workspace
        badchans = input('Bad channels? [f.i. [1 5 9] or leave blank]:');
        if badchans
            % remove bad channels
            EEG = pop_select(EEG,'nochannel',badchans);
            EEG.badchan = badchans;
        end
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step1.set',eegfile(1:end-4)),'filepath','.');
        
        % reconstitute regions matrix [start,stops] from EEG boundary events
        [starts,stops] = reconstitute_regions(EEG.event);
        %[starts,stops] = combineregions([starts' stops']); % remove overlapping regions
        fid = fopen(sprintf('%s_step1_rejected.txt',eegfile(1:end-4)),'wt');
        for i = 1:numel(starts)
            fprintf(fid,'%i\t%i\n',starts(i),stops(i));
        end
        fclose(fid);
        close_down

    case 2
        disp('Step 2: ICA of pruned data');
        EEG = pop_loadset(sprintf('%s_step1.set',eegfile(1:end-4)));
        % run ICA on pruned data
        EEG = pop_runica(EEG,'extended',1,'interupt','off');
        EEG.setname='Step2result';
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step2.set',eegfile(1:end-4)),'filepath','.');
        %pop_expica(EEG,'weights','/Users/marc/Documents/MATLAB/marc_bus1_ICA.txt');
    
    case 3
        disp('Step 3: Prune ICA activations');
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        EEG = pop_loadset(sprintf('%s_step2.set',eegfile(1:end-4)));
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
        EEG.setname='Step3result';
        % rename old boundary events
        for i = 1:length(EEG.event)
            if strcmp(EEG.event(i).type,'boundary')
                EEG.event(i).type = 'boundary_old';
            end
        end
        pop_eegplot(EEG,0,1,1);
        assignin('base','EEG',EEG); % move EEG variable to base workspace because pop_eegplot is not blocking
        uiwait;
        EEG = evalin('base','EEG'); % recover modified EEG variable from base workspace
        
        % reconstitute regions matrix [start,stops] from EEG boundary events
        [starts,stops] = reconstitute_regions(EEG.event);
        %[starts,stops] = combineregions([starts' stops']); % remove overlapping regions
        fid = fopen(sprintf('%s_step3_rejected.txt',eegfile(1:end-4)),'wt');
        for i = 1:numel(starts)
            fprintf(fid,'%i\t%i\n',starts(i),stops(i));
        end
        fclose(fid);
        
        % revert renaming of old boundary events
        for i = 1:length(EEG.event)
            if strcmp(EEG.event(i).type,'boundary_old')
                EEG.event(i).type = 'boundary';
            end
        end
        
        % save pruned data set
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step3.set',eegfile(1:end-4)),'filepath','.');
        close_down;
        
    case 4
        disp('Step 4: Merge data sets from same session');
        [FileName,PathName] = uigetfile('.set','Select SET files to merge:','multiselect','on');
        nfiles = length(FileName);
        if ~iscell(FileName); error('No files selected!'); end
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        for i = 1:nfiles
            EEG = pop_loadset([PathName FileName{i}]);
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,0);
        end
        EEG = pop_mergeset(ALLEEG,[1:nfiles],0);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step4.set',[eegfile(1:end-4)]),'filepath','.');
        
        % save list of files for later use
        fid = fopen(sprintf('%s_step4_filelist.txt',eegfile(1:end-4)),'wt');
        for i = 1:1:nfiles
            fprintf(fid,'%s ',FileName{i});
        end
        fclose(fid);
        close_down;
        
    case 5
        disp('Step 5: ICA on merged data for clustering');
        EEG = pop_loadset(sprintf('%s_step4.set',eegfile(1:end-4)));
        EEG = pop_runica(EEG,'extended',1,'interupt','off');
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step5.set',eegfile(1:end-4)),'filepath','.');
        
    case 6
        disp('Step 6: Restart from raw data and reject datapoints (ensures consistent rejections in EEG and audio data)');
        % go through all files from one session
        fid = fopen(sprintf('%s_step4_filelist.txt',eegfile(1:end-4)),'r');
        line = fgetl(fid);
        fclose(fid);
        files = regexp(line,' ','split');
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        for current_file = 1:length(files)-1
            file = files{current_file}(1:end-10);
            eegfile = sprintf('%s.bdf',file);
            fprintf('loading %s...\n',eegfile);
            % preprocess again
            EEG = pop_biosig(eegfile,'channels',[1:3,8:22]);
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,0);
            EEG.chanlocs = readlocs('chanlocs.ced');
            EEG.data(16,:) = EEG.data(16,:)/2;
            EEG = pop_reref(EEG,16);
            EEG = pop_select(EEG,'nochannel',11);
            EEG2 = pop_loadset(sprintf('%s_step1.set',file));
            if isfield(EEG2,'badchan') % remove same bad channels as in step 1
                EEG.badchan = EEG2.badchan;
                EEG = pop_select(EEG,'nochannel',EEG2.badchan);
            end
            clear EEG2;
            EEG = pop_eegfiltnew(EEG,[],1,8250,true,[],0);
            EEG = pop_eegfiltnew(EEG,[],40,330,0,[],0);
            EEG.setname='Step6result';
            len1 = EEG.pnts;
            fprintf('Initial length of EEG recording: %i. \n',len1);
            % load deletions
            rej_1 = load(sprintf('%s_step1_rejected.txt',file));
            rej_2 = load(sprintf('%s_step3_rejected.txt',file));
            % construct vector of non-rejected sample points
            remaining_samples = 1:size(EEG.data,2);
            deleted = [];
            for i = 1:size(rej_1,1);
                deleted = [deleted rej_1(i,1):rej_1(i,2)];
            end
            del1 = length(deleted);
            fprintf('Deleting %i samples in first round, ',del1);
            remaining_samples(deleted) = [];
            deleted = [];
            for i = 1:size(rej_2,1);
                deleted = [deleted rej_2(i,1):rej_2(i,2)];
            end
            del2 = length(deleted);
            fprintf('and %i samples in second round.\n',del2);
            remaining_samples(deleted) = [];
            % remove sample points from EEG.data
            EEG.data = EEG.data(:,remaining_samples);
            EEG.pnts = size(EEG.data,2);
            EEG = eeg_checkset(EEG);
            len2 = EEG.pnts;
            fprintf('Length after deletions: %i. \n',len2);
            assert(len1-del1-del2 == len2); % make sure that the deletion add up, otherwise abort
            [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,current_file); % write changes back into ALLEEG (necessary? check!)
            EEG = pop_saveset(EEG,'filename',sprintf('%s_step6a.set',file),'filepath','.');
        end
        % Merge files
        EEG = pop_mergeset(ALLEEG,1:length(files)-1,0);
        % Add IC decomposition from step 5
        EEG2 = pop_loadset(sprintf('%s_step5.set',files{1}(1:end-10)));
        EEG.icaweights = EEG2.icaweights;
        EEG.icawinv = EEG2.icawinv;
        EEG.icasphere = EEG2.icasphere;
        clear EEG2;
        EEG = eeg_checkset(EEG);
        fprintf('Final length of merged EEG recordings: %i. \n',EEG.pnts);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step6.set',files{1}(1:end-10)),'filepath','.');
        close_down;
        
    case 7
        disp('Step 7: Select ICs');
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        EEG = pop_loadset(sprintf('%s_step6.set',eegfile(1:end-4)));
        if isempty(EEG.chanlocs(1).X); % no coordinates loaded
            EEG = pop_chanedit(EEG,'lookup','chanlocs.ced');
        end
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG,EEG,0);
        ncomp = size(EEG.icaweights,1);
        % plot activations
        pop_eegplot(EEG,0,0,0);
        % plot IC topo at 10 Hz
        %dur = (size(EEG.icaact,2)/EEG.srate)*1000-2; % data length in ms
        %figure; pop_spectopo(EEG,0,[0 dur],'EEG','freq',[10],'plotchan',0,'percent',100,'icacomps',[1:ncomp],'nicamaps',ncomp,'freqrange',[1 30],'electrodes','off');
        % topoplot all ICs
        pop_selectcomps(EEG,[1:ncomp]);
        assignin('base','EEG',EEG); % move EEG variable to base workspace because pop_eegplot is not blocking
        assignin('base','ALLEEG',ALLEEG); %
        uiwait(gcf);
        EEG = evalin('base','EEG'); % recover modified EEG variable from base workspace
        ALLEEG = evalin('base','ALLEEG');
        % rejected ICs now in EEG.reject.gcompreject
        EEG = eeg_checkset(EEG,'ica');
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step6.set',eegfile(1:end-4)),'filepath','.');
        % save accepted IC activations
        good_ICs = find(EEG.reject.gcompreject==0);
        icaact = EEG.icaact(good_ICs,:);
        save(sprintf('%s_step6_icaact',eegfile(1:end-4)),'icaact','good_ICs');
        close_down;
        
    case 8
        disp('Step 7: Save IC spectrogram features.');
        disp('If component clustering was performed, then this is done on the group ICs,');
        disp('otherwise on the individual ICs.');
        if exist(sprintf('.%s%s_groupICs.txt',filesep,eegfile(1:end-4)),'file') % group ICs available
            groupICs = load(sprintf('.%s%s_groupICs.txt',filesep,eegfile(1:end-4))); % load the list
            EEG = pop_loadset(sprintf('%s_step6.set',eegfile(1:end-4))); % load EEG data to extract ICs listed in group IC file
            icaact = EEG.icaact(groupICs,:);
        else
            EEG = pop_loadset(sprintf('%s_step6.set',eegfile(1:end-4))); % load EEG data to extract ICs listed in group IC file
            load(sprintf('%s_step6_icaact',eegfile(1:end-4))); % now we have icaact in memory
        end
        % compute the spectrogram of each IC, taking boundaries into account
        % extract boundaries
        bounds = .5; % this will be sample 1 after removing offset and adding 1 (see below)
        for i = 1:length(EEG.event)
            if strcmp(EEG.event(i).type,'boundary')
                bounds = [bounds EEG.event(i).latency];
            end
        end
        bounds = bounds - 0.5;   % remove offset added by eegrej.m to move boundary between samples
        bounds = bounds + 1; % eegrej subtracts 1 at the start
        % compute spectrogram for each data segment (not crossing boundaries to avoid edge artifacts) 
        icaact_specgram = []; % container for resulting spectrograms, ic x time x freq
        for i = 1:numel(bounds)-1
            win = bounds([i-1 i]);
            spec = dBspectrogram(icaact(:,win(1):win(2)));
            icaact_specgram = cat(2,icaact_specgram,spec); % concatenate sprecgrams time axis
        end
        
end


function close_down
close all;
clear all;


function [starts,stops] = reconstitute_regions(events)
% list of deletions (start and stops in samples)
del = [];
n = 1;
for i = 1:length(events)
    if strcmp(events(i).type,'boundary') && ~isnan(events(i).duration)
        starts(n) = events(i).latency;
        durs(n) = events(i).duration;
        n = n + 1;
    end
end
starts = starts - 0.5;   % remove offset added by eegrej.m to move boundary between samples
for i = 1:numel(starts) % add cut-out durations back to star times
    for idx = i+1:numel(starts)
        starts(idx) = starts(idx) + durs(i);
    end
end
stops = starts + durs; % calculate stop times from durations and starts
starts = starts + 1; % eegrej subtracts 1 at the start


function [starts,stops] = combineregions(regions)
% copied from eeglab's eeg_eegrej.m
regions = sortrows(sort(regions,2)); % Sorting regions
allreg = [regions(:,1)' regions(:,2)'; ones(1,numel(regions(:,1))) -ones(1,numel(regions(:,2)')) ].';
allreg = sortrows(allreg,1); % Sort all start and stop points (column 1),

mboundary = cumsum(allreg(:,2)); % Rationale: regions will start always with 1 and close with 0, since starts=1 end=-1
indx = 0; count = 1;

while indx ~= length(allreg) 
    starts(count) = allreg(indx+1,1);
    [tmp,I]= min(abs(mboundary(indx+1:end)));
    stops(count) = allreg(I + indx,1);
    indx = indx + I ;
    count = count+1;
end

function spec = dBspectrogram(sig,w)
% calculate spectra of all rows (ICs) in sig 
% sig  ... ic x time
% w    ... window (such as hamming(1024))
% spec ... ic x time x freq

hN = (N/2)+1; % size of positive spectrum, including sample 0
hM1 = floor((length(w)+1)/2); % half analysis window size by rounding
hM2 = floor(length(w)/2); % half analysis window size by floor
fftbuffer = zeros(N,1); % initialize buffer for FFT
w = w/sum(w); % normalize analysis window
tol = 1e-10; % tolearance for setting bins to zero; improves phase angle calculation

%%% TEST THIS! then add support for 2D
for i = 1:size(sig,1)
    xw = x.*w; % window the input sound
    fftbuffer(1:hM1) = xw(hM2:end); % zero-phase window in fftbuffer
    fftbuffer(end-hM2:end) = xw(1:hM2);
    X = fft(fftbuffer); % compute FFT
    X = X(1:hN); % only use positive frequencies
    absX = abs(X); % compute absolute value of positive side
    absX(absX<eps) = eps; % if zeros add epsilon to handle log
    mX = 20 * log10(absX); % magnitude spectrum of positive frequencies in dB
    X(real(abs(real(X))) < tol) = 0.0; % for phase calculation set to 0 the small values
    X(imag(abs(imag(X))) < tol) = 0.0; % for phase calculation set to 0 the small values         
    pX = unwrap(angle(X)); % unwrapped phase spectrum of positive frequencies
end

%EEG = pop_subcomp(EEG,[1 2 3],0);
% concatenate all recordings, then do ICA on filtered and pruned data
% (1-40/50 Hz), (then perhaps reject again on ICA time course), then apply
% weights to pruned raw data (although we don#t actually do that because we
% will predict from the ICA time courses directly).
% Also, calculate contralateralty and ITD/ILD relationship!
% Cross-validation to see whether that is predictive
