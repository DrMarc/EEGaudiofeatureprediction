function preprocess_features(eegfile);

disp('Remove rejected data from features, concatenate runs and z-score');
disp('The result should be aligned features and ICA ready for machine learning');

fid = fopen(sprintf('%s_step4_filelist.txt',eegfile(1:end-4)),'r');
line = fgetl(fid);
fclose(fid);
files = regexp(line,' ','split');

for current_file = 1:length(files)-1
    file = files{current_file}(1:end-10);
    load(sprintf('%s-features.mat',file));
    % load deletions
    rej_1 = load(sprintf('%s_step1_rejected.txt',file));
    rej_2 = load(sprintf('%s_step3_rejected.txt',file));
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
    % save as new feature vector, all in one large matrix (samples x features)
    % also save a cell array of feature names and frequencies corresponding
    % to the rows of the feature vector
    
    % Before we do that, a hack to removen an incosistency in the way the
    % 1D features are saved by combineEEGaudifeatures.m
    % This should be taken out and fixed in combineEEGaudifeatures.m
    features.Elhilali.mapped2EEG.saliency = features.Elhilali.mapped2EEG.saliency';
    features.Elhilali.mapped2EEG.envelope = features.Elhilali.mapped2EEG.envelope';
    features.Elhilali.mapped2EEG.pitch = features.Elhilali.mapped2EEG.pitch';
    features.Elhilali.mapped2EEG.specgram = features.Elhilali.mapped2EEG.specgram';
    features.Elhilali.mapped2EEG.bw = features.Elhilali.mapped2EEG.bw';
    features.Elhilali.mapped2EEG.rate = features.Elhilali.mapped2EEG.rate';
    features.audio.mapped2EEG.spectral_centroid = features.audio.mapped2EEG.spectral_centroid';
    features.audio.mapped2EEG.spectral_brightness = features.audio.mapped2EEG.spectral_brightness';
    features.audio.mapped2EEG.spectral_flux = features.audio.mapped2EEG.spectral_flux';
    % %%%
    if current_file == 1
        [featurevec,featurenames] = lcf_getfeatures(features,remaining_samples);
    else
        [tmpvec,tmpnames] = lcf_getfeatures(features,remaining_samples);
        if isequal(tmpnames,featurenames)
            % concatenate
            featurevec = [featurevec;tmpvec];
        else
            featurenames, tmpnames
            error('Features not identical across files!');
        end
    end
end
% z-score all features
for i = 1:size(featurevec,2)
    nmean = nanmean(featurevec(:,i));
    nstd  = nanstd(featurevec(:,i));
    featurevec(:,i) = (featurevec(:,i)-nmean)/nstd;
end
save(sprintf('%s_featurevec.mat',eegfile(1:end-4)),'featurevec','featurenames');
fprintf('Saved %s_featurevec.mat\n',eegfile(1:end-4));


function [featurevec,featurenames] = lcf_getfeatures(features,remaining_samples)
highlevel_fields = fieldnames(features); % should be {'Elhilali','Kayser','audio'}
featurenames = {};
pos = 1;
for f = 1:length(highlevel_fields);
    fields = eval(sprintf('fieldnames(features.%s.mapped2EEG)',highlevel_fields{f}));
    for i = 1:length(fields);
        if eval(sprintf('isstruct(features.%s.mapped2EEG.%s)',highlevel_fields{f},fields{i})) % necessary because of 'norm' field (struct, but no feature)
            continue;
        else
            len = eval(sprintf('min(size(features.%s.mapped2EEG.%s))',highlevel_fields{f},fields{i}));
            if len > 1
                eval(sprintf('featurevec(:,pos:pos+len-1) = features.%s.mapped2EEG.%s(remaining_samples,:);',highlevel_fields{f},fields{i}));
                for k = 0:len-1
                    featurenames{pos+k} = sprintf('%s %i',fields{i},eval(sprintf('round(features.%s.freqs(k+1))',highlevel_fields{f})));
                end 
            else
                eval(sprintf('featurevec(:,pos) = features.%s.mapped2EEG.%s(remaining_samples,:);',highlevel_fields{f},fields{i}));
                featurenames{pos} = fields{i};
            end
            pos = pos + len;
        end
    end % end of loop through all audio fields in one run 
end
        


