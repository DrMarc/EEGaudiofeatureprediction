function EEGanalysis(eegfile)

% determine step
step = 1;
while exist(sprintf('.%s%s_step%i.set',filesep,eegfile(1:end-4),step),'file');
    step = step + 1;
end

switch step
    case 1
        global EEG
        disp('Step 1: Manual pruning of continuous data');
        EEG = pop_biosig(eegfile,'channels',[1:3,8:22]);
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
        [start,stops] = reconstitute_regions(EEG.event);
        fid = fopen(sprintf('%s_step1_rejected.txt',eegfile(1:end-4)),'wt');
        fprintf(fid,'%i\t%i\n',start,stops);
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
        EEG = pop_loadset(sprintf('%s_step2.set',eegfile(1:end-4)));
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
        [start,stops] = reconstitute_regions(EEG.event);
        fid = fopen(sprintf('%s_step3_rejected.txt',eegfile(1:end-4)),'wt');
        fprintf(fid,'%i\t%i\n',start,stops);
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
        disp('Step 6: Select ICs');
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        EEG = pop_loadset(sprintf('%s_step5.set',eegfile(1:end-4)));
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
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG,'filename',sprintf('%s_step6.set',eegfile(1:end-4)),'filepath','.');
        % save accepted IC activations
        good_ICs = find(EEG.reject.gcompreject==0);
        icaact = EEG.icaact(good_ICs,:);
        save(sprintf('%s_step6_icaact',eegfile(1:end-4)),'icaact','good_ICs');
        close_down;
        

end

function close_down
close all;
clear all;

function newregions = combineregions(regions)
% copied from eeglab's eeg_eegrej.m
regions = sortrows(sort(regions,2)); % Sorting regions
allreg = [ regions(:,1)' regions(:,2)'; ones(1,numel(regions(:,1))) -ones(1,numel(regions(:,2)')) ].';
allreg = sortrows(allreg,1); % Sort all start and stop points (column 1),

mboundary = cumsum(allreg(:,2)); % Rationale: regions will start always with 1 and close with 0, since starts=1 end=-1
indx = 0; count = 1;

while indx ~= length(allreg) 
    newregions(count,1) = allreg(indx+1,1);
    [tmp,I]= min(abs(mboundary(indx+1:end)));
    newregions(count,2) = allreg(I + indx,1);
    indx = indx + I ;
    count = count+1;
end

% Verbose
if size(regions,1) ~= size(newregions,1)
    warning('Overlapping regions detected and fixed.');
end

function [start,stops] = reconstitute_regions(events)
% list of deletions (start and stops in samples)
del = [];
n = 1;
for i = 1:length(events)
    if strcmp(events(i).type,'boundary') && ~isnan(events(i).duration)
        start(n) = events(i).latency;
        durs(n) = events(i).duration;
        n = n + 1;
    end
end
start = start - 0.5;   % remove offset added by eegrej.m to move boundary between samples
for i = 1:numel(start) % add cut-out durations back to star times
    for idx = i+1:numel(start)
        start(idx) = start(idx) + durs(i);
    end
end
stops = start + durs; % calculate stop times from durations and starts
start = start + 1; % eegrej subtracts 1 at the start

%EEG = pop_subcomp(EEG,[1 2 3],0);
% concatenate all recordings, then do ICA on filtered and pruned data
% (1-40/50 Hz), (then perhaps reject again on ICA time course), then apply
% weights to pruned raw data (although we don#t actually do that because we
% will predict from the ICA time courses directly).
% Also, calculate contralateralty and ITD/ILD relationship!
% Cross-validation to see whether that is predictive
