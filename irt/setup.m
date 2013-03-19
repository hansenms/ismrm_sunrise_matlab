% setup.m
% run this file to set up matlab path etc.
% you may need to modify this depending on how you installed the toolbox
% so this should be considered simply a "guide" not a robust script.

if ~exist('irtdir', 'var')
	disp('The variable "irtdir" is not set, so trying default, assuming')
	disp('that you launched matlab from the irt install directory.')
	disp('You may need to edit setup.m or adjust your path otherwise.')

%	irtdir = pwd; % the default is to assume launch from irt directory

	% default is to look for directory where this setup.m is installed!
	irtdir = which('setup'); % find setup.m
	[irtdir dummy] = fileparts(irtdir);
	clear dummy

	disp(['Assuming you installed irt in directory "' irtdir '".'])

%	irtdir = '~fessler/l/src/matlab/alg/'; % where you install this package
%	irtdir = '~fessler/l/web/irt/'; % where you install this package
end

if ~exist(irtdir, 'dir')
	disp(sprintf('The directory "%s" does not exist', irtdir))
	error(sprintf('you need to edit %s to change default path', mfilename))
end

if irtdir(end) ~= filesep % make sure there is a '/' at end of directory
	irtdir = [irtdir filesep];
end

addpath([irtdir 'align']);		% image registration
addpath([irtdir 'align/mex']);		% image registration mex files
addpath([irtdir 'blob']);		% blob (KB) basis
addpath([irtdir 'ct']);			% x-ray CT (polyenergetic) recon
addpath([irtdir 'data']);		% example data
addpath([irtdir 'emission']);		% emission image reconstruction
addpath([irtdir 'example']);		% example applications
addpath([irtdir 'fbp']);		% FBP (filtered backprojection) code
addpath([irtdir 'general']);		% generic image reconstruction
addpath([irtdir 'graph']);		% graphics routines
addpath([irtdir 'mri']);		% MRI reconstruction
addpath([irtdir 'mri-rf/yip-spsp']);	% MRI RF pulse design
%addpath([irtdir 'mri/recon']);		% MRI reconstruction - old
addpath([irtdir 'nufft']);		% nonuniform FFT (for a fast projector)
addpath([irtdir 'nufft/table']);	% mex files for NUFFT
addpath([irtdir 'penalty']);		% regularizing penalty functions
addpath([irtdir 'systems']);		% system "matrices"
addpath([irtdir 'systems/tests']);	% tests of systems
addpath([irtdir 'transmission']);	% transmission image reconstruction
addpath([irtdir 'utilities']);		% various utility functions
addpath([irtdir 'wls']);		% weighted least-squares (WLS) estimates

tmp = [irtdir 'um']; % extra files for um users only
if exist(tmp, 'dir'), addpath(tmp), end
tmp = [irtdir 'um/lustig-poisson-disk'];
if exist(tmp, 'dir'), addpath(tmp), end

%if isequal(version, '3.6.3') % octave
if ir_is_octave
	addpath([irtdir 'octave']); % extra stuff for octave only!
elseif isempty(which('dbstack')) % for freemat only!
	addpath([irtdir 'freemat']); % extra stuff for freemat only!
end

% Set up path to mex files, possibly depending on matlab version.
% Fortunately it seems that the v6 ones will run on v7 too.
% Unfortunately, it seems that sometimes compiling on one version
% e.g., 7.3 and running on earlier 7.0 won't work.
% if you have this problem, then comment out the mex path line(s) below.
% Much of the toolbox will work without mex, just slower.
%if str2num(version('-release')) <= 13
%	addpath([irtdir 'mex/v6'])
%else
	addpath([irtdir 'mex/v7']);
%end

% check to see if path setup worked by looking for im() routine.
if strcmp([irtdir 'graph' filesep 'im.m'], which('im'))
	disp('Path setup for irt appears to have succeeded.')
else
	disp('Path setup for irt may have failed.')
end
