% dye_diffusion_simulation.m v1.0.0
% This script was developed by Agata Witkowska for the study
% Witkowska A. & Jahn R., Biophysical Journal (2017)
% "Rapid SNARE-mediated fusion of liposomes and chromaffin granules with giant unilamellar vesicles"
% http://dx.doi.org/10.1016/j.bpj.2017.03.010

% This is a simulation of labelled lipid diffusion in the membrane upon fusion of a vesicle containing such labelled lipid.
% Additional details are provided in accompanying .pdf file.
% It runs with Matlab and GNU Octave
% Required GNU Octave packages: statistics.

clear all;
close all;

% Test for Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave == true
  % Needed for Octave
  pkg load statistics;
  graphics_toolkit gnuplot;
end

% Parameters set
N = 985; % number of dye particles in the vesicle
D = 2.83; % mean value obtained from FRAP on GUVs with TR-PE in um^2/s (Fig. 1 in the main paper)
dt = 0.001; % time step given in s
radius = 0.1; % radius of the circle
binSize = 0.005; % bin size for histograms
runMax = 50; % number of runs for statistics
tPoints = 50; % number of time points
snapshots = [5,10,20,50]; % tPoints to plot

% Check input snapshots
if isempty(snapshots) || (size(snapshots,2) == 1 && any(snapshots == 0))
  error('At least one snapshot !=0 is needed');
end
if ~any(snapshots==0)
  snapshots = [snapshots 0]; % Snapshot at 0 is needed
end
snapshots = sort(snapshots); % Snapshots sorting

name = @(snapshot) genvarname(strjoin({'R',num2str(snapshot)},''));

for snapshot = snapshots
  data.(name(snapshot)) = [];
end;

% Generation of a waitbar
h = waitbar (0, '', 'createcancelbtn', 'setappdata(gcbf,''canceling'',1)');

dx = sqrt(dt * 2 * D); % approximation of particle displacement in 1D according to Einstein [Einstein, A. (1905). Über die von der molekularkinetischen Theorie der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen. Annalen der Physik 322, 549–560.]

for run = 1:runMax
  % Initial random distribution of N particles in a circular spot coming from fusion of a liposome
  angle = 2 * pi * rand(1, N);
  rr = radius * sqrt(rand(1, N));
  r = [rr .* cos(angle); rr .* sin(angle)];
  
  % snapshot at 0
  data.(name(0)) = [data.(name(0)), r];

  for t = 1:tPoints
    r = r + dx * randn(2, N);
    if any(snapshots == t)
      data.(name(t)) = [data.(name(t)), r];
    end
  end

  % Check for Cancel button press
  if getappdata(h, 'canceling')
    break
  end
  waitbar (run / runMax, h, ['Run ' num2str(run) ' of ' num2str(runMax)]);

end
delete(h)

% Plot generation
for snapshot = snapshots
  figure
  if snapshot == 0
    maxVal = radius;
  else
    maxVal = max(max(abs(data.(name(snapshot)))));
  end
  H0x = hist3(data.(name(snapshot))', [round(maxVal/binSize) round(maxVal/binSize)]);
  if snapshot == 0
    H00 = H0x;
    H00c = [zeros(size(H00, 1), 1) H00 zeros(size(H00, 1), 1)];
	  interval = size(H00c, 2);
    n1 = H00';
  else
	  interval = size(H0x, 1);
    n1 = H0x / (N * run)';
  end
  n1(size(H0x, 1) + 1, size(H0x, 2) + 1) = 0;
  xb = linspace(min(data.(name(snapshot))(1, :)), max(data.(name(snapshot))(1, :)), size(H0x, 1) + 1);
  yb = linspace(min(data.(name(snapshot))(2, :)), max(data.(name(snapshot))(2, :)), size(H0x, 1) + 1);
  h = pcolor(xb, yb, n1);
  if snapshot ~= 0
    maximum = max(max(H00 / length(data.(name(0)))));
	  caxis([0, maximum]); % adjusts colormap display to the histogram at t=0
  end
  colormap(hot);
  set(h, 'edgecolor', 'none');
  grid off;
  axis square;
  axis([- 1.5, 1.5, - 1.5, 1.5]);
  set(gca, 'color', [10 / 255 0 0]);

  print([num2str(snapshot) 'ms'], '-dsvg');

  figure
  BinNumber = round(round(maxVal / binSize) / 2);
  if snapshot == 0
    yline = H00c(BinNumber, :) / (N * run);
  else
    yline = H0x(BinNumber, :) / (N * run);
  end
  xline = linspace(min(data.(name(snapshot))(1, :)), max(data.(name(snapshot))(1, :)), interval);
  plot(xline, yline);
  if snapshot <= snapshots(2) % normalize line plots y-display
    axisLimit = max(H0x(BinNumber, :) / (N * run));
  end
  axis square;
  axis([- 1.5 1.5 0 axisLimit]);

  print([num2str(snapshot) 'msL'], '-dsvg');
end

TotalTimeInMs = t * dt * 1000 % simulation time in ms