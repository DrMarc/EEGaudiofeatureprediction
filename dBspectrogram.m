function X = dBspectrogram(x)

fs = 500; % sampling rate in Hz
dbdown = 100; % cutoff for faint components (-100 dB)
M = 100; % window length 200 ms @ 500 Hz FS
Mhalf = M/2; % half of the window
w = 0.54 - 0.46 * cos(2*pi*(0:M-1)'/M); % causal Hamming window
w = w/(M*0.54); % scale for unity overlap-add with overlap of M-1
%nfft = 2^(nextpow2(M)); % will be 128 - not using it, i.e. no spectral interpolation
freqs = (0:Mhalf-1)/Mhalf * fs/2; % frequency vector
freqs = freqs(2:11); % only use useful frequencies below 50 Hz (data is filtered), exclude DC
nhop = 1; % advance one sample per frame
tol = 1e-10; % tolearance for setting bins to zero; improves phase angle calculation
x = x(:); % make sure it's a column
if length(x)<M % zero-pad to fill a window:
  x = [x;zeros(M-length(x),1)];
end;
nx = length(x);
nframes = nx;
%X = zeros(Mhalf,nframes); % allocate output spectrogram
X = zeros(10,nframes); % allocate output spectrogram, but only save useful freqs below 50 Hz, excluding DC
%PX = zeros(Mhalf,nframes); % allocate output phase
xframe = zeros(M,1); 
xoff = -Mhalf; % input time offset = half a frame
for m = 1:nframes
  if xoff<0 % the window will start initially before the signal, only get partial singal then (will be zero-padded)
    xframe(1:xoff+M) = x(1:xoff+M); % partial input data frame
  else
    if xoff+M > nx
      xframe = [x(xoff+1:nx);zeros(xoff+M-nx,1)]; % at the end of the signal, windows will be padded with zeroes
    else
      xframe = x(xoff+1:xoff+M); % input data frame
    end
  end
  xw = w .* xframe; % Apply window
  xwzp = [xw(Mhalf+1:M);xw(1:Mhalf)]; % moving the zero point
  frame_spec = fft(xwzp);
  %frame_spec = frame_spec(1:Mhalf); % only use positive frequencies
  frame_spec = frame_spec(2:11); % actually, only use useful frequencies below 50 Hz (data is filtered), exclude DC
  absX = abs(frame_spec)*2; % compute absolute value of positive side and multiply by 2 to account for trowing out the negative frequencies
  absX(absX<eps) = eps; % if zeros add epsilon to handle log
  X(:,m) = 20 * log10(absX); % magnitude spectrum of positive frequencies in dB
  %frame_spec(abs(real(frame_spec)) < tol) = 0.0; % for phase calculation set to 0 the small values
  %frame_spec(abs(imag(frame_spec)) < tol) = 0.0; % for phase calculation set to 0 the small values         
  %PX(:,m) = unwrap(angle(frame_spec)); % unwrapped phase spectrum of positive frequencies
  xoff = xoff + nhop; % advance input offset by hop size
end

if (nargout==0)
  fs = 500;
  t = (0:nframes-1)*nhop/fs;
  Xmax = max(max(X));
  % Clip lower limit to -dbdown dB so nulls don't dominate:
  clipvals = [Xmax-dbdown,Xmax];
  imagesc(t,freqs,X,clipvals);
  % grid;
  axis('xy');
  colormap(jet);
  xlabel('Time (sec)');
  ylabel('Freq (Hz)');
end