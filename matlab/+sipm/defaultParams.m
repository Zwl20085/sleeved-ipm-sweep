function s = defaultParams(topology, varargin)
%SIPM.DEFAULTPARAMS Default parameters of a topology as a MATLAB struct.
%
%   s = sipm.defaultParams('sleeved_bridged')
%   s = sipm.defaultParams('sleeved_bridgeless', 'sleeve_thickness', 1.5)
%
% Topologies: 'sleeved_bridged', 'sleeved_bridgeless'.
% The Python side validates the values (an invalid set raises an error).
% Units: mm, deg, GPa, kg/m^3, rpm, A, degC, 1/K (see docs/parameters.md).
% The fields are read from the Python package (sipm params --json), so new parameters
% appear here without changes to this file.

if nargin < 1, topology = 'sleeved_bridged'; end
sets = '';
for k = 1:2:numel(varargin)
    sets = [sets, sprintf(' %s=%s', varargin{k}, toText(varargin{k+1}))]; %#ok<AGROW>
end
if ~isempty(sets), sets = [' --set', sets]; end
out = sipm.cli(sprintf('params --json --topology %s%s', topology, sets));
lines = strsplit(strtrim(out), newline);
s = jsondecode(lines{end});
end

function t = toText(v)
if ischar(v) || isstring(v)
    t = ['"\"', char(v), '\""'];
else
    t = sprintf('%.17g', v);
end
end
