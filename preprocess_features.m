function preprocess_features(eegfile);

disp('Remove rejected data from features, concatenate runs and z-score');
disp('The result should be aligned features and ICA ready for machine learning');

fid = fopen(sprintf('%s_step4_filelist.txt',eegfile(1:end-4)),'r');
line = fgetl(fid);
fclose(fid);
files = regexp(line,' ','split');

for i = 1:length(files)-1
    eegfile = files{i}(1:end-10);
    load(sprintf('%s-features.mat',eegfile));
    % load deletions
    rej_1 = load(sprintf('%s_step1_rejected.txt',eegfile));
    rej_2 = load(sprintf('%s_step3_rejected.txt',eegfile));
    % construct vector of non-rejected sample points
    remaining_samples = 1:length(features.audio.mapped2EEG.spectral_centroid);
    for i = 1:size(rej_1,1);
        deleted = [1:rej_1(i,2)]+rej_1(i,1)-1;
        remaining_samples = setdiff(remaining_samples,deleted);
    end
    for i = 1:size(rej_2,1);
        deleted = [1:rej_2(i,2)]+rej_2(i,1)-1;
        remaining_samples = setdiff(remaining_samples,deleted);
    end
    % delete those samples from each entry in features.X.mapped2EEG and
    % save as new feature vectors
    featurevec(:,1:16) = features.audio.mapped2EEG.ITD(remaining_samples,:);
end

featurevec(:,17:) = features.audio.mapped2EEG.ILD(remaining_samples,:);
features_trimmed.audio.onsets = features.audio.mapped2EEG.onsets(remaining_samples,:);
features_trimmed.audio.offsets = features.audio.mapped2EEG.offsets(remaining_samples,:);
features_trimmed.audio.spectral_centroid = features.audio.mapped2EEG.spectral_centroid(remaining_samples);
features_trimmed.audio.spectral_brightness = features.audio.mapped2EEG.spectral_brightness(remaining_samples);
features_trimmed.audio.spectral_flux = features.audio.mapped2EEG.spectral_flux(remaining_samples);

features_trimmed.Kayser.saliency = features.audio.mapped2EEG.ITD(remaining_samples,:);
        


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

