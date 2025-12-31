function setup()
% Add package folders to path for current session.
% Usage: run `setup` once per MATLAB session or add to your project startup.

  here = fileparts(mfilename('fullpath'));
  addpath(here, fullfile(here, '+sgdelta'));
  if exist('tests','dir'), addpath(fullfile(here,'tests')); end
  fprintf('[sgdelta] Path added. Use functions as sgdelta.construct_delta_allpass(...)\n');
end
