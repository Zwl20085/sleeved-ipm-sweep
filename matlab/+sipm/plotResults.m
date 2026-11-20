function files = plotResults(resultsCsv, physics, outDir, opts)
%SIPM.PLOTRESULTS Standard figure set (PNG) for a results CSV written by a sweep.
%
%   files = sipm.plotResults(fullfile(tempdir, 'sipm_runs', 'demo', 'results.csv'), 'mech', 'figs', ...
%                            'SigmaAllow', 400, 'InterferenceFraction', 0.2)
%
% physics: 'em', 'mech' or 'acloss'. For 'mech', SigmaAllow (MPa) enables the
% interpolated stress map with its feasible region. Returns the files written.

arguments
    resultsCsv (1,:) char
    physics (1,:) char {mustBeMember(physics, {'em', 'mech', 'acloss'})}
    outDir (1,:) char
    opts.SigmaAllow double = []
    opts.InterferenceFraction (1,1) double = 0.2
end
args = sprintf('plot "%s" --physics %s --out "%s"', resultsCsv, physics, outDir);
if ~isempty(opts.SigmaAllow)
    args = sprintf('%s --sigma-allow %.17g --interference-fraction %.17g', args, opts.SigmaAllow, ...
        opts.InterferenceFraction);
end
out = sipm.cli(args);
tok = regexp(out, 'figure: ([^\r\n]+)', 'tokens');
files = cellfun(@(c) strtrim(c{1}), tok, 'UniformOutput', false)';
end
