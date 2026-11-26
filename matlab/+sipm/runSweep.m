function T = runSweep(name, physics, design, opts)
%SIPM.RUNSWEEP Run (or resume) a sweep and return the results as a MATLAB table.
%
%   T = sipm.runSweep(name, physics, design)
%   T = sipm.runSweep(name, physics, design, 'Topology', 'sleeved_bridged', ...
%                     'WorkDir', fullfile(tempdir, 'sipm_runs'), 'Base', struct('mech_rpm', 20000))
%
%   name     folder name of the sweep under WorkDir
%   physics  'em' (average torque at gamma, current), 'mech' (stresses) or 'acloss'
%   design   MATLAB table, one row per point, variable names = parameter names
%            (see sipm.defaultParams); sipm.grid builds a full factorial
%
% Name-value options:
%   Topology     default topology for rows without a 'topology' column
%   Base         struct of overrides applied to every row
%   WorkDir      run directory (default: SIPM_WORKDIR or <tempdir>/sipm_runs)
%   StopAfter    stop after this many newly solved points (resume by calling again)
%   RetryFailed  re-run points that failed before (default false)
%   AedtVersion  e.g. '2024.1' (default: newest installed)
%   MapdlCores   cores for MAPDL (default 2)
%   MinFreeRamGb      wait before each point until this much RAM is free (default 4)
%   RamWaitTimeoutS   then stop the sweep cleanly; call again to resume (default 1800)
%   MaxSolverRamGb    terminate a solver above this RSS, point recorded as failed (default: off)
%   AedtRestartEvery  release and restart AEDT every N EM points (default 10, 0 = never)
%
% Each point is cached on disk as soon as it finishes; calling the same sweep
% again skips cached points. Failed and invalid points are returned with their
% message in the 'message' column instead of stopping the sweep. The sweep is
% described by a TOML file written next to the results, so it can be re-run with
% the command line:  sipm sweep <WorkDir>/<name>/matlab_sweep.toml

arguments
    name (1,:) char
    physics (1,:) char {mustBeMember(physics, {'em', 'mech', 'acloss'})}
    design table
    opts.Topology (1,:) char = 'sleeved_bridged'
    opts.Base struct = struct()
    opts.WorkDir (1,:) char = ''
    opts.StopAfter double = []
    opts.RetryFailed (1,1) logical = false
    opts.AedtVersion (1,:) char = ''
    opts.MapdlCores (1,1) double = 2
    opts.MinFreeRamGb (1,1) double = 4
    opts.RamWaitTimeoutS (1,1) double = 1800
    opts.MaxSolverRamGb double = []
    opts.AedtRestartEvery (1,1) double = 10
end

workdir = opts.WorkDir;
if isempty(workdir), workdir = getenv('SIPM_WORKDIR'); end
if isempty(workdir), workdir = fullfile(tempdir, 'sipm_runs'); end
sdir = fullfile(workdir, name);
if ~isfolder(sdir), mkdir(sdir); end

L = {'[sweep]', kv('name', name), kv('physics', physics), kv('topology', opts.Topology), ...
     kv('workdir', strrep(workdir, '\', '/')), '', '[ansys]', kv('mapdl_cores', opts.MapdlCores)};
if ~isempty(opts.AedtVersion), L{end+1} = kv('aedt_version', opts.AedtVersion); end
L{end+1} = kv('min_free_ram_gb', opts.MinFreeRamGb);
L{end+1} = kv('ram_wait_timeout_s', opts.RamWaitTimeoutS);
L{end+1} = kv('aedt_restart_every', opts.AedtRestartEvery);
if ~isempty(opts.MaxSolverRamGb), L{end+1} = kv('max_solver_ram_gb', opts.MaxSolverRamGb); end
L{end+1} = '';
L{end+1} = '[base]';
f = fieldnames(opts.Base);
for j = 1:numel(f), L{end+1} = kv(f{j}, opts.Base.(f{j})); end %#ok<AGROW>
vars = design.Properties.VariableNames;
for i = 1:height(design)
    L{end+1} = ''; %#ok<AGROW>
    L{end+1} = '[[points]]'; %#ok<AGROW>
    for j = 1:numel(vars)
        v = design{i, j};
        if iscell(v), v = v{1}; end
        L{end+1} = kv(vars{j}, v); %#ok<AGROW>
    end
end
cfgFile = fullfile(sdir, 'matlab_sweep.toml');
fid = fopen(cfgFile, 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);

args = sprintf('sweep "%s" --no-plots', cfgFile);
if ~isempty(opts.StopAfter), args = sprintf('%s --stop-after %d', args, opts.StopAfter); end
if opts.RetryFailed, args = [args, ' --retry-failed']; end
out = sipm.cli(args, true);
if contains(out, 'STOPPED:')
    warning('sipm:stopped', 'Sweep stopped by the memory gate; call sipm.runSweep again to resume.');
end
T = readtable(fullfile(sdir, 'results.csv'), 'TextType', 'string', 'Delimiter', ',');
end

function s = kv(k, v)
if ischar(v) || isstring(v)
    s = sprintf('%s = "%s"', k, char(v));
elseif islogical(v)
    s = sprintf('%s = %s', k, lower(string(v)));
elseif v == round(v) && abs(v) < 1e15
    s = sprintf('%s = %d', k, v);
else
    s = sprintf('%s = %.17g', k, v);
end
end
