function compare_marker_files(audiomarkerfile,EEGmarkerfile)
% read both files
timestampsaudio = load(audiomarkerfile);
timestampsEEG = [];
n = 1;
fid = fopen(EEGmarkerfile,'r');
while ~feof(fid)
    line = fgetl(fid);
    while line
        [tmp,line] = strtok(line);
    end
    timestampsEEG(n) = str2num(tmp);
    n = n+1;
end
fclose(fid);

% % translate to common time base (secs from start of recording)
% timestampsAudio_sec = timestampsAudio./FS_audio;
% timestampsEEG_sec = timestampsEEG./FS_EEG;
% disp('Please check timestamps!');
% keyboard; % if they don't line up, do sth like >>timestampsEEG = timestampsEEG([1 2 4 6 7 8]);

