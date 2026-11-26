%EXAMPLE_TINY_SWEEP Tiny mechanical sweep driven entirely from MATLAB.
%
% Runs 2 sleeve thicknesses x 3 interferences for the sleeved bridge-retained rotor
% at 20000 rpm through MAPDL, returns a MATLAB table and plots it. Needs a licensed
% MAPDL; about 10-20 s per point. Set SIPM_PYTHON to the python.exe of the
% environment where sleeved-ipm-sweep is installed (or pass it to sipm.setup), and
% optionally SIPM_WORKDIR for the run directory.
%
% Run from the repository root:   matlab -batch "run('matlab/examples/example_tiny_sweep.m')"

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));             % the folder that contains +sipm
sipm.setup();

p = sipm.defaultParams('sleeved_bridged');
fprintf('lamination OD %.2f mm for a %.1f mm sleeve (air gap %.2f mm)\n', ...
    p.stator_id - 2*p.airgap - 2*p.sleeve_thickness, p.sleeve_thickness, p.airgap);

design = sipm.grid('sleeve_thickness', [1.0 2.0], 'interference', [0.0 0.1 0.2]);
T = sipm.runSweep('matlab_tiny_mech', 'mech', design, 'Topology', 'sleeved_bridged');
disp(T(:, {'sleeve_thickness', 'interference', 'status', 'svm_rotor_max', 'sh_sleeve_max', 'wall_s'}));

ok = T(T.status == "ok", :);
f = figure('Visible', 'off');
tiledlayout(1, 2);
nexttile; hold on; grid on;
for t = unique(ok.sleeve_thickness)'
    r = ok(ok.sleeve_thickness == t, :);
    plot(r.interference, r.svm_rotor_max, '-o', 'DisplayName', sprintf('%.1f mm sleeve', t));
end
xlabel('Radial interference (mm)'); ylabel('Max von Mises, lamination outside magnet contacts (MPa)'); legend('Location', 'best');
nexttile; hold on; grid on;
for t = unique(ok.sleeve_thickness)'
    r = ok(ok.sleeve_thickness == t, :);
    plot(r.interference, r.sh_sleeve_max, '-s', 'DisplayName', sprintf('%.1f mm sleeve', t));
end
xlabel('Radial interference (mm)'); ylabel('Max sleeve hoop stress (MPa)'); legend('Location', 'best');
workdir = getenv('SIPM_WORKDIR');
if isempty(workdir), workdir = fullfile(tempdir, 'sipm_runs'); end
png = fullfile(workdir, 'matlab_tiny_mech', 'figures', 'matlab_stress_lines.png');
if ~isfolder(fileparts(png)), mkdir(fileparts(png)); end
exportgraphics(f, png);
fprintf('MATLAB figure: %s\n', png);

% the same results through the Python plotting functions
files = sipm.plotResults(fullfile(workdir, 'matlab_tiny_mech', 'results.csv'), 'mech', ...
    fullfile(workdir, 'matlab_tiny_mech', 'figures'), 'SigmaAllow', 400);  % example yield, use your material's
fprintf('Python figure: %s\n', files{:});
