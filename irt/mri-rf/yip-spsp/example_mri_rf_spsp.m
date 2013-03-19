% example_mri_rf_spsp
% This script is an example of spectral-spatial pulse design for
% through-plane phase precompensatory slice selection for T2*-weighted
% functional MR, based on our publication:
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectral-spatial RF pulse design script, based on "Spectral-spatial RF
% pulse design for through-plane phase precompensatory slice selection for
% T2*-weighted functional MRI", Chun-yu Yip et al, Magnetic Resonance in
% Medicine, 2009. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by Chun-yu Yip, University of Michigan, Ann Arbor, 4/1/2009
%
% Running this script requires that the entire image reconstruction toolbox
% by Professor Jeffrey A. Fessler be installed first: 
%
% http://www.eecs.umich.edu/~fessler/code/index.html
%
% You can conveniently download the whole package at
%
% http://www.eecs.umich.edu/~fessler/irt/fessler.tgz
%
% and use the unix command 'tar' to retrieve the files and folders.
% Use the "setup.m" command in the toolbox to set up paths first.


% First load pulse design parameters into matlab structures kp, rfp, iop.
if ~isvar('kp'), printm 'Loading parameters...'
	[kp iop rfp] = example_spsp_param1();
end

iop.waveformfilespath = './out/'; % jf: put output in a subdirectory
if ~exist(iop.waveformfilespath, 'dir')
	fail('need to create output directory "%s"', iop.waveformfilespath)
end

% Design z-gradient waveform, based on parameters in kp.
if ~isvar('gz'), printm 'Designing z gradient waveform...'
	[kp gz kz kf t] = compute_gz_spsp(kp);
end

if 1 % display gradient waveform
	clf, pl = @(i) subplot(230+i);
	pl(1)
	plot(t*1000, gz), axis tight
	xlabel('time (msec)')
	ylabel('g/cm')
	title('Gz waveform')
	grid
end

% Design complex-valued RF waveform iteratively using conjugate gradient
if ~isvar('b'), printm 'Computing RF pulse waveform...'
	[b d f z mm] = compute_rf_spsp_mgh(kp, rfp, gz, kz, kf);

	% Write computed waveforms to files for simulation and/or scanner.
%	iop.writetofile_sim = false; % uncomment to save file i/o
%	iop.writetofile_scanner = false;
	write2files_spsp(kp, rfp, iop, gz, b);
end

ir_mri_rf_spsp_plot(b, gz, d, f, z, mm, t) % display results
prompt

% Perform Bloch simulation in SPSP space
if ~isvar('mresult'), printm 'Performing Bloch simulation...'
	mresult = dosim7_spsp(kp, rfp, iop);
end
