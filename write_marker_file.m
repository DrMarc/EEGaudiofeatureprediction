function write_marker_file(file)

type = file(end-2:end);

switch type
    case 'wav'
        markerfile = [file(1:end-4) '-wavmarkers.txt'];
        if exist(['.' filesep markerfile],'file')
            if ~strcmpi(input('Overwrite exisiting file [y/n]?','s'),'y')
                error('Exiting without changes')
            end
        end
        % read the audio data
        [~,~,~,opt_ck] = my_wavread(file); % using the extended wavread that can load CUE points
        timestampsAudio = opt_ck.cue_sampleoffset;
        fid = fopen(markerfile,'wt');
        fprintf(fid,'%i\n',timestampsAudio);
        fclose(fid);
    
    case 'bdf'
        markerfile = [file(1:end-4) '-markers.txt'];
        if exist(['.' filesep markerfile],'file')
            if ~strcmpi(input('Overwrite exisiting file [y/n]?','s'),'y')
                error('Exiting without changes')
            end
        end
        EEG = pop_biosig(file,'channels',[]);
        Nmarkers = length(EEG.event);
        for i = 1:Nmarkers
            timestampsEEG(i) = EEG.event(i).latency; % in samples from start
        end
        fid = fopen(markerfile,'wt');
        fprintf(fid,'1, Marker 1, %i\n',timestampsEEG); % "1, Marker 1, " is added to that the format is the same as in the native marker files 
        fclose(fid);
    otherwise
        error('Allowed file types are wav and bdf.');
end
