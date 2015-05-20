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

% show a formatted list for comparison
lcf_printlist(timestampsaudio,timestampsEEG);

% match lists
if length(timestampsaudio)<=length(timestampsEEG) % prune longer EEG list
    timestampsEEG = lcf_matchlists(timestampsaudio,timestampsEEG);
else % prune longer audio list
    timestampsaudio = lcf_matchlists(timestampsEEG,timestampsaudio);
end

timestampsaudio_bkp = timestampsaudio;
timestampsEEG_bkp = timestampsEEG;
response = 'n';
while response ~= 'y'
    disp('Matched points:');
    lcf_printlist(timestampsaudio,timestampsEEG);
    response = input('Write this to file? [y/n]','s');
    if response == 'y'
        lcf_writelists(timestampsaudio,timestampsEEG,audiomarkerfile,EEGmarkerfile);
    else
        disp('Select matching points (ex. timepointsEEG = timepointsEEG([1,2,4:8]))');
        disp('or restore original: timepointsEEG = timepointsEEG_bkp');
        disp('or type "dbquit" to exit without changes.');
        keyboard;
    end
end


function time2 = lcf_matchlists(time1,time2)
% normalize both lists to [0..1]
t1 = time1-time1(1);
t1 = t1/t1(end);
t2 = time2-time2(1);
t2 = t2/t2(end);
% go through the first list (presumed to be shorter) and find matching timepoints (nearest neighbor),
% keep a list of matched points
matched_points = [];
for point = 1:length(t1)
    [~,matched_points(point)] = min(abs(t2-t1(point)));
end
% discard unmatched points from longer list
time2 = time2(matched_points);


function lcf_printlist(timestampsaudio,timestampsEEG)
% translate to common time base (secs from start of recording)
timestampsaudio_sec = timestampsaudio./48000;
timestampsEEG_sec = timestampsEEG./500;
timestampsaudio_sec = timestampsaudio_sec - timestampsaudio_sec(1);
timestampsEEG_sec = timestampsEEG_sec - timestampsEEG_sec(1);
% show list of time points
fprintf(1,'Nr\taudio\tEEG\n');
n = length(timestampsaudio);
m = length(timestampsEEG);
for i = 1:max(n,m)
    fprintf(1,'%i\t',i);
    if i<=n
        fprintf(1,'%.5g\t',timestampsaudio_sec(i));
    else
        fprintf(1,'\t');
    end
    if i<=m
        fprintf(1,'%.5g\t',timestampsEEG_sec(i));
    else
        fprintf(1,'\t');
    end
    fprintf(1,'\n');
end


function lcf_writelists(timestampsaudio,timestampsEEG,audiomarkerfile,EEGmarkerfile)
% write audio
fid = fopen(audiomarkerfile,'wt');
fprintf(fid,'%i\n',timestampsaudio);
fclose(fid);
% write EEG
fid = fopen(EEGmarkerfile,'wt');
fprintf(fid,'1, Marker 1, %i\n',timestampsEEG); % "1, Marker 1, " is added to that the format is the same as in the native marker files 
fclose(fid);