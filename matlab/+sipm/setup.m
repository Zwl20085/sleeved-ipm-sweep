function exe = setup(pythonExe)
%SIPM.SETUP Tell the MATLAB front end which Python has sleeved-ipm-sweep installed.
%
%   sipm.setup('C:\path\to\env\python.exe')
%   sipm.setup()            % uses the SIPM_PYTHON environment variable
%
% The MATLAB functions run the Python package as a separate process (its command
% line interface), so Python, Ansys and MATLAB never share DLLs or OpenMP settings.
% The interpreter is stored in the SIPM_PYTHON environment variable of this
% MATLAB session.

if nargin < 1 || isempty(pythonExe)
    pythonExe = getenv('SIPM_PYTHON');
end
if isempty(pythonExe)
    error('sipm:setup', 'Pass the python.exe of your environment or set SIPM_PYTHON.');
end
if ~isfile(pythonExe)
    error('sipm:setup', 'Python executable not found: %s', pythonExe);
end
setenv('SIPM_PYTHON', pythonExe);
[st, out] = system(sprintf('"%s" -s -c "import sipm; print(sipm.__version__)"', pythonExe));
if st ~= 0
    error('sipm:setup', ['Cannot import the Python package "sipm" with %s.\n' ...
        'Install it with:  python -m pip install -e <repo>[aedt]\n%s'], pythonExe, out);
end
fprintf('sleeved-ipm-sweep %s via %s\n', strtrim(out), pythonExe);
exe = pythonExe;
end
