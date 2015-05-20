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
timestampsaudio_sec = timestampsaudio./48000;
timestampsEEG_sec = timestampsEEG./500;
timestampsaudio_sec = timestampsaudio_sec - timestampsaudio_sec(1);
timestampsEEG_sec = timestampsEEG_sec - timestampsEEG_sec(1);
plot(timestampsaudio_sec,'ro');
hold on;
plot(timestampsEEG_sec,'b.');

% normalize both lists to [0..1]
tmp1 = timestampsaudio-timestampsaudio(1);
tmp1 = tmp1/tmp1(end);
tmp2 = timestampsEEG-timestampsEEG(1);
tmp2 = tmp2/tmp2(end);
% go through the longer list and find a timepoint that does not have a
% matching point in the other list

