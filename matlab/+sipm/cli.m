function out = cli(args, echo)
%SIPM.CLI Run the sipm command line interface with the configured Python.
%
%   out = sipm.cli('params --json --topology sleeved_bridged')
%
% Errors if the command returns a non-zero status. With echo = true the output is
% streamed to the command window while the command runs.

if nargin < 2, echo = false; end
py = getenv('SIPM_PYTHON');
if isempty(py)
    error('sipm:cli', 'Call sipm.setup(pythonExe) first (or set SIPM_PYTHON).');
end
cmd = sprintf('"%s" -s -m sipm.cli %s', py, args);   % -s: ignore the per-user site folder
if echo
    [st, out] = system(cmd, '-echo');
else
    [st, out] = system(cmd);
end
if st ~= 0 && ~(echo && any(st == [2 3]))   % 2 = some points failed, 3 = stopped by the RAM gate
    error('sipm:cli', 'sipm %s failed (status %d):\n%s', args, st, out);
end
end
