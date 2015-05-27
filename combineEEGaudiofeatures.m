function features = combineEEGaudiofeatures(EEGfile,audiofile)
% iterates through an audio recording, extracts features and timestamps,
% and resamples features to align with EEG data

% read the audio data
[audio,FS_audio,~,opt_ck] = my_wavread(audiofile); % using the extended wavread that can load CUE points
if exist(['.' filesep audiofile(1:end-4) '-wavmarkers.txt'],'file') % read markers from text file
    timestampsAudio = load(['.' filesep audiofile(1:end-4) '-wavmarkers.txt']);
else
    timestampsAudio = opt_ck.cue_sampleoffset;
end
decimation_factor = 6; % decimate a 48k recording to 8k
tmp_r = decimate(audio(:,1),decimation_factor);
tmp_l = decimate(audio(:,2),decimation_factor);
audio = single([tmp_r,tmp_l]);
clear tmp_r tmp_l;
FS_audio = FS_audio/decimation_factor;
% How do we keep the timestamps intact? Just divide the sample values by the decimation factor.
timestampsAudio = timestampsAudio/decimation_factor;

% Now get the EEG timestamps
% the goal is to resample features so that they align with the EEG and add
% them as channels to the EEG
[timestampsEEG,FS_EEG,NsamplesEEG] = lcf_getEEGtimestamps(EEGfile);

% feature vectors and time codes ready to resample and combine with EEG
features.Elhilali = lcf_extractaudiofeatures_saliencyElhilali(audio(:,1),FS_audio);
features.Kayser = lcf_extractaudiofeatures_saliencyKayser(audio(:,1),FS_audio);
features.audio = lcf_extractaudiofeatures_spatial(audio,FS_audio);

% convert timestamps from audio to each feature time course
features.audio.timestamps = timestampsAudio/(FS_audio/features.audio.FS); % same sample rate
features.Kayser.timestamps = timestampsAudio/(FS_audio/features.Kayser.FS); % 
features.Elhilali.timestamps = timestampsAudio/(FS_audio/features.Elhilali.FS); % same sample rate

% Now get the EEG timestamps
% the goal is to resample features so that they align with the EEG and add
% them as channels to the EEG
[timestampsEEG,FS_EEG,NsamplesEEG] = lcf_getEEGtimestamps(EEGfile);

% treat the timestamp values as x and y coordinates, then a regression line
% gives the optimal mapping from feature samples to EEG samples 
p = polyfit(timestampsEEG,timestampsAudio,1);
% calculate mean deviation as sanity check
fitted_timestampsAudio = polyval(p,timestampsEEG);
standart_deviation = sqrt(sum(((fitted_timestampsAudio - timestampsAudio)/FS_audio).^2)) % in sec
figure;
plot(timestampsEEG/FS_EEG,timestampsAudio/FS_audio,'bo',timestampsEEG/FS_EEG,fitted_timestampsAudio/FS_audio,'r-');
xlabel('timestamps EEG');
ylabel('timestamps Audio');

% resample features to fit with EEG
features.audio.p = polyfit(timestampsEEG,features.audio.timestamps,1);
features.Kayser.p = polyfit(timestampsEEG,features.Kayser.timestamps,1);
features.Elhilali.p = polyfit(timestampsEEG,features.Elhilali.timestamps,1);
% resampling the feature vectors to fit with the EEG, starting with audio features
mappedSamplesFeatues2EEG = polyval(features.audio.p,1:NsamplesEEG);
% 2D features require interp2, resample only along the time axis
Y = 1:size(features.audio.ITD,1);
X = (1:size(features.audio.ITD,2))';
features.audio.mapped2EEG.ITD = interp2(features.audio.ITD,X,mappedSamplesFeatues2EEG,'cubic',NaN);
features.audio.mapped2EEG.ILD = interp2(features.audio.ILD,X,mappedSamplesFeatues2EEG,'cubic',NaN);
features.audio.mapped2EEG.onsets = interp2(features.audio.onsets,X,mappedSamplesFeatues2EEG,'cubic',NaN);
features.audio.mapped2EEG.offsets = interp2(features.audio.offsets,X,mappedSamplesFeatues2EEG,'cubic',NaN);
features.audio.mapped2EEG.spectral_centroid = interp1(Y,features.audio.spectral_centroid,mappedSamplesFeatues2EEG,'cubic',NaN);
features.audio.mapped2EEG.spectral_brightness = interp1(Y,features.audio.spectral_brightness,mappedSamplesFeatues2EEG,'cubic',NaN);
features.audio.mapped2EEG.spectral_flux = interp1(Y,features.audio.spectral_flux,mappedSamplesFeatues2EEG,'cubic',NaN);
% audio features are smaple accurate

% same for Kayser saliency
mappedSamplesFeatues2EEG = polyval(features.Kayser.p,1:NsamplesEEG);
% 2D features require interp2, resample only along the time axis
X = (1:size(features.Kayser.saliency',2))';
features.Kayser.mapped2EEG.saliency = interp2(features.Kayser.saliency',X,mappedSamplesFeatues2EEG,'cubic',NaN);
% Kayser features appear to be shifted backward in time by 0.05 sec!

% same for Elhilali features (all 1D)
mappedSamplesFeatues2EEG = polyval(features.Elhilali.p,1:NsamplesEEG);
Y = 1:size(features.Elhilali.saliency',1);
features.Elhilali.mapped2EEG.saliency = interp1(Y,features.Elhilali.saliency',mappedSamplesFeatues2EEG,'cubic',NaN);
features.Elhilali.mapped2EEG.envelope = interp1(Y,features.Elhilali.envelope',mappedSamplesFeatues2EEG,'cubic',NaN);
features.Elhilali.mapped2EEG.pitch = interp1(Y,features.Elhilali.pitch',mappedSamplesFeatues2EEG,'cubic',NaN);
features.Elhilali.mapped2EEG.specgram = interp1(Y,features.Elhilali.specgram',mappedSamplesFeatues2EEG,'cubic',NaN);
features.Elhilali.mapped2EEG.bw = interp1(Y,features.Elhilali.bw',mappedSamplesFeatues2EEG,'cubic',NaN);
features.Elhilali.mapped2EEG.rate = interp1(Y,features.Elhilali.rate',mappedSamplesFeatues2EEG,'cubic',NaN);
% Elhilali features appear to be shifted forward in time by 2 samples!
% Elhilali saliency has 3 peaks with 2-sample-or-so diff before onset of
% salient tone! 

% calculate and save scaling factors to normalize standard dev of each
% feature to 1; these are stand-ins or the normalized features that will be calculated when needed
features.audio.mapped2EEG.norm.ITD = [-nanmean(features.audio.mapped2EEG.ITD(:)) 1/nanstd(features.audio.mapped2EEG.ITD(:))];
features.audio.mapped2EEG.norm.ILD = [-nanmean(features.audio.mapped2EEG.ILD(:)) 1/nanstd(features.audio.mapped2EEG.ILD(:))];
features.audio.mapped2EEG.norm.onsets = [-nanmean(features.audio.mapped2EEG.onsets(:)) 1/nanstd(features.audio.mapped2EEG.onsets(:))];
features.audio.mapped2EEG.norm.offsets = [-nanmean(features.audio.mapped2EEG.offsets(:)) 1/nanstd(features.audio.mapped2EEG.offsets(:))];
features.audio.mapped2EEG.norm.spectral_centroid = [-nanmean(features.audio.mapped2EEG.spectral_centroid) 1/nanstd(features.audio.mapped2EEG.spectral_centroid)];
features.audio.mapped2EEG.norm.spectral_brightness = [-nanmean(features.audio.mapped2EEG.spectral_brightness) 1/nanstd(features.audio.mapped2EEG.spectral_brightness)];
features.audio.mapped2EEG.norm.spectral_flux = [-nanmean(features.audio.mapped2EEG.spectral_flux) 1/nanstd(features.audio.mapped2EEG.spectral_flux)];
features.Kayser.mapped2EEG.norm.saliency = [-nanmean(features.Kayser.mapped2EEG.saliency(:)) 1/nanstd(features.Kayser.mapped2EEG.saliency(:))];
features.Elhilali.mapped2EEG.norm.saliency = [-nanmean(features.Elhilali.mapped2EEG.saliency) 1/nanstd(features.Elhilali.mapped2EEG.saliency)];
features.Elhilali.mapped2EEG.norm.envelope = [-nanmean(features.Elhilali.mapped2EEG.envelope) 1/nanstd(features.Elhilali.mapped2EEG.envelope)];
features.Elhilali.mapped2EEG.norm.pitch = [-nanmean(features.Elhilali.mapped2EEG.pitch) 1/nanstd(features.Elhilali.mapped2EEG.pitch)];
features.Elhilali.mapped2EEG.norm.specgram = [-nanmean(features.Elhilali.mapped2EEG.specgram) 1/nanstd(features.Elhilali.mapped2EEG.specgram)];
features.Elhilali.mapped2EEG.norm.bw = [-nanmean(features.Elhilali.mapped2EEG.bw) 1/nanstd(features.Elhilali.mapped2EEG.bw)];
features.Elhilali.mapped2EEG.norm.rate = [-nanmean(features.Elhilali.mapped2EEG.rate) 1/nanstd(features.Elhilali.mapped2EEG.rate)];


%%% local functions %%%

function features = lcf_extractaudiofeatures_spatial(audio,FS)
disp('Computing audio features...');
startAuditoryFrontEnd
dataObj = dataObject(audio,FS,floor(length(audio)/FS)); % dangerous - buffering full audio because the chunking doesn't seem to work
% Request interaural time differences (ITDs)
requests = {'itd','ild','onset_map','offset_map','spectral_features'};
% Parameters of the auditory filterbank processor
fb_type       = 'gammatone';
fb_lowFreqHz  = 100;
fb_highFreqHz = 4000;
fb_nChannels  = 16;  
% Parameters of innerhaircell processor
ihc_method    = 'dau';
% Parameters of crosscorrelation processor
cc_wSizeSec  = 0.02;
cc_hSizeSec  = 0.01;
cc_wname     = 'hann';
% Parameters of ratemap processor
rm_wSizeSec  = 0.02;
rm_hSizeSec  = 0.01;
rm_decaySec  = 8E-3;
rm_wname     = 'hann';
% Parameters of transient detector (same parameters for onsets & offests)
trm_on_minValuedB    = -80;
trm_on_minStrengthdB = 3;
trm_on_minSpread     = 5;
trm_on_fuseWithinSec = 30E-3;
% Summary of parameters 
par = genParStruct('fb_type',fb_type,'fb_lowFreqHz',fb_lowFreqHz,...
                   'fb_highFreqHz',fb_highFreqHz,'fb_nChannels',fb_nChannels,...
                   'ihc_method',ihc_method,'cc_wSizeSec',cc_wSizeSec,...
                   'cc_hSizeSec',cc_hSizeSec,'cc_wname',cc_wname,'ac_wSizeSec',rm_wSizeSec,...
                   'ac_hSizeSec',rm_hSizeSec,'rm_decaySec',rm_decaySec,...
                   'ac_wname',rm_wname,'trm_minValuedB',trm_on_minValuedB,...
                   'trm_minStrengthdB',trm_on_minStrengthdB,'trm_minSpread',trm_on_minSpread,...
                   'trm_fuseWithinSec',trm_on_fuseWithinSec); 
managerObj = manager(dataObj,requests,par);
managerObj.processSignal(); % Request processing
% postprocess features
freqs = dataObj.gammatone{1}.cfHz;
ITDs = dataObj.itd{1}.Data(:);
features.ITD = ITDs;
ILDs = dataObj.ild{1}.Data(:);
ILDs(ILDs<-10) = -10; % clip at 10dB (20dB range) to limi the influence of artifacts
ILDs(ILDs>10) = 10;
features.ILD = ILDs;
onsets  = dataObj.onset_map{1}.Data(:);
features.onsets = onsets;
offsets = dataObj.offset_map{1}.Data(:);
features.offsets = offsets;
tmp = dataObj.spectral_features{1}.Data(:);
features.spectral_centroid  = tmp(:,1);
features.spectral_brightness = tmp(:,3);
features.spectral_flux      = tmp(:,13);
features.FS = dataObj.itd{1}.FsHz;
features.freqs = freqs;
features.t = (0:length(features.ITD)-1)/features.FS;


function features = lcf_extractaudiofeatures_saliencyKayser(audio,FS_audio)
disp('Computing Kayser saliency...');
spec_window = 800;
spec_overlap = 780;
nx = length(audio); % size of signal
FS_out = (FS_audio/((spec_window-spec_overlap)))/2; % output FS depends on specgram window and overlap (divide by 2 because specgram is subsampled at factor 2)
[spec,freqs,tw] = specgram(audio,1024,FS_audio,spec_window,spec_overlap); % compute spectrogram
freqs = freqs(1:2:end)'+eps; % saliency map is stored at half the resolution of the specgram
nfreqs = length(freqs); % number of freqs in each window
tw = tw(1:2:end); % saliency map is stored at half the time resolution of the specgram
spec = log(abs(spec)); % make intensity map
tmp = Saliency_map(spec,4); % compute saliency map
features.saliency = tmp.eo + tmp.esi + tmp.epi;
features.FS = FS_out;
features.freqs = freqs; % saliency map is stored at half the resolution of the specgram


function features = lcf_extractaudiofeatures_saliencyElhilali(audio,FS_audio)
%load testData.mat; %loads data, fs_orig, W
window_len = 10; % sec
cfg = [];
cfg.fs_resample = FS_audio;
%cfg.W = W;
cfg.W = eye(5); % perhaps safer to use this
cfg.downsample_factor = 16; % this determines the output fs
disp('Computing Elhilali saliency...');
% divide up a signal into windows
nx = length(audio); % size of signal
nw = FS_audio * window_len; % size of window in samples
overlap = FS_audio * .5; % 1/2s of overlap
nwindows = floor(nx/(nw-overlap));
L_frm = round(1 * 2^(4+log2(FS_audio/16000))); %length of spectrogram windows
FS_out = FS_audio/L_frm/cfg.downsample_factor; %fs adjusted for all the downsampling for saliency estimate
ST = 1/FS_out;
t = 0:ST:window_len-ST;
ntw = length(t); % number of time values in each output window
overlap_out = round(FS_out * 0.25); % 15s of overlap at output time resolution
nx_out = (nwindows*(ntw-overlap_out))+overlap_out; % total length of audio over window length - overlap, then add one overlap back to account for the end = total output timesteps
pos = 1;
pos_out = 1;
sal = zeros(1,nx_out);
streams = zeros(5,nx_out); % Envelope, Pitch, Spectrogram, Bandwidth, Rate
fprintf(1,'This will take ~%i min.',round((nx/428)/60)); %
n = 1 % loop counter
while (pos+nw-1 <= nx) % while enough signal left
    y = audio(pos:pos+nw-1); % make window y
    [tmp_sal,tmp_streams] = get_sal(y,FS_audio,cfg);
    tmp_streams = cell2mat(tmp_streams); % convert cell array to vector;
    tmp_streams = reshape(tmp_streams,ntw,5)'; % reshape 
    % add raised cosine ramps
    tmp_sal = lcfCos2Gate(tmp_sal,overlap_out);
    tmp_streams = lcfCos2Gate(tmp_streams,overlap_out);
    sal(pos_out:pos_out+ntw-1) = sal(pos_out:pos_out+ntw-1) + tmp_sal';
    streams(:,pos_out:pos_out+ntw-1) = streams(:,pos_out:pos_out+ntw-1) + tmp_streams;
    pos = pos + nw - overlap; % next window
    pos_out = pos_out+ntw-overlap_out; % next output window
    n = n+1 % loop counter
end
features.saliency = sal;
features.envelope = streams(1,:);
features.pitch = streams(1,:);
features.specgram = streams(1,:);
features.bw = streams(1,:);
features.rate = streams(1,:);
features.FS = FS_out;

function [timestampsEEG,FS_EEG,NsamplesEEG] = lcf_getEEGtimestamps(EEGfile)
EEG = pop_biosig(EEGfile,'channels',[]);
FS_EEG = EEG.srate;
NsamplesEEG = EEG.pnts;
if exist(['.' filesep EEGfile(1:end-4) '-markers.txt'],'file') % read markers from text file
    disp('Reading markers from text file.')
    timestampsEEG = [];
    n = 1;
    fid = fopen(['.' filesep EEGfile(1:end-4) '-markers.txt'],'r');
    while ~feof(fid)
        line = fgetl(fid);
        while line
            [tmp,line] = strtok(line);
        end
        timestampsEEG(n) = str2num(tmp);
        n = n+1;
    end
    fclose(fid);
else
    warning('No text file. Reading markers from BDF file.')
    Nmarkers = length(EEG.event);
    for i = 1:Nmarkers
        timestampsEEG(i) = EEG.event(i).latency; % in samples (FS_EEG) from start
    end
end

function out = flatten(in,freqs,aweighting)
% Comput the mean across fregs with or without A-weighting.
% in is assumed to be a 2d matrix of time x freq
if aweighting
    Ra = @(f) (12200^2*f.^4)./( (f.^2*20.6^2) .* sqrt((f.^2+107.7^2).*(f.^2+737.9^2)) .* (f.^2 + 12200^2) ); % A-weights as anonymous function
    w = Ra(freqs); % weight for each frequency
    w = w/sum(w); % normalize weights to unity
    w = repmat(w,size(in,1),1); % expand to the size of in
    in = in .* w; % apply weights
end
out = mean(in,2);

function out = lcfCos2Gate(in,GPts)
% adds ramps of length GPts (in samples) to a wave form
if min(size(in))>1 % 2D array; gating along 2. dimension!
    SPts = size(in,2);
    Nfreqs = size(in,1);
    if SPts<2*GPts
        error('sound is shorter than 2 * GPts!')
    end
    env = 0.5 - 0.5 * cos(0:pi/GPts:pi); % raised cosine window
    env = repmat(env,Nfreqs,1);
    out = in .* [env ones(Nfreqs,SPts-size(env,2)*2) fliplr(env)];
else % 1D (vector)
    SPts = length(in);
    if SPts<2*GPts
        error('sound is shorter than 2 * GPts!')
    end
    env = 0.5 - 0.5 * cos(0:pi/GPts:pi)'; % raised cosine window
    out = in .* [env;ones(length(in)-length(env)*2,1);flipud(env)];
end
