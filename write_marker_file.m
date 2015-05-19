function write_marker_file(audiofile)

markerfile = [audiofile(1:end-4) '-wavmarkers.txt'];
% read the audio data
[~,~,~,opt_ck] = my_wavread(audiofile); % using the extended wavread that can load CUE points
timestampsAudio = opt_ck.cue_sampleoffset;
if exist(['.' filesep markerfile],'file')
    if ~strcmpi(input('Overwrite exisiting file [y/n]?','s'),'y')
        error('Exiting without changes')
    end
end
fid = fopen(markerfile,'wt');
fprintf(fid,'%i\n',timestampsAudio);
fclose(fid);