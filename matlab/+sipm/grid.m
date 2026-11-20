function T = grid(varargin)
%SIPM.GRID Full-factorial design table.
%
%   T = sipm.grid('sleeve_thickness', [1 2], 'interference', [0 0.1 0.2])
%
% Returns a MATLAB table with one row per combination, ready for sipm.runSweep.

names = varargin(1:2:end);
vals = varargin(2:2:end);
n = numel(names);
grids = cell(1, n);
[grids{:}] = ndgrid(vals{:});
T = table();
for k = 1:n
    T.(names{k}) = grids{k}(:);
end
T = sortrows(T, names);
end
