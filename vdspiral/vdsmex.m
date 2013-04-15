%
%	function [k,g,s,time] = vdsmex(N,FOV,res,gmax,smax,T,ngmax)
%
%
%	VARIABLE DENSITY SPIRAL GENERATION:
%	----------------------------------
%	This function uses a mex file to quickly generate a
%	variable-density spiral.  See vds.m for more details.
%
%	INPUT:
%		N	Number of interleaves.
%		FOV	Field-of-View* 	 (cm)
%		res	Resolution.	 (mm)
%		gmax	max gradient amplitude (G/cm)	[3.9]
%		smax 	maximum slew rate (G/cm/s)	[14500]
%		T 	sampling period (s)		[0.000004]
%		ngmax   max # points along gradient.    [10000000]
%	OUTPUT:
%		k	k-space points (Ng x 2)
%		g	gradient waveform (Ng x 2)
%		s	slew waveform (Ng x 2)
%		time	time waveform (Ng x 1)
%
%	type vdsmex.m for more info...

%	Note that FOV can be a vector, in which case it is
%	interpreted as the polynomial coefficients of the FOV
%	in terms of k-space radius, ie
%
%	field of view =    sum[i=1:length(FOV)]  FOV(i)*(kr/krmax)^(i-1) 
%
%	The spiral ends either when the desired resolution (res)
%	is reached, or when the maximum number of gradient points
%	(ngmax) is reached.

% =============== CVS Log Messages ==========================
%	$Log: vdsmex.m,v $
%	Revision 1.3  2004/09/08 20:57:13  brian
%	--Fixed comment on FOV expression.
%	--Added ngmax as a parameter to vdsmex, since vds_mex.c already supports it.
%	--Fixed bug where ngmax is now multiplied by oversamp so that
%	  ngmax specifies the OUTPUT gradient length.
%	
%	Revision 1.2  2003/05/29 23:04:19  brian
%	Modified to support the
%	max points OR max k-space radius option.
%	
%	Revision 1.1  2003/04/23 21:46:45  brian
%	Added mex version of vds.m
%	
%	Revision 1.3  2002/11/18 05:36:02  brian
%	Rounds lengths to a multiple of 4 to avoid
%	frame size issues later on.
%	
%	Revision 1.2  2002/11/18 05:32:19  brian
%	minor edits
%	
%	Revision 1.1  2002/03/28 01:03:20  bah
%	Added to CVS
%	
%
% ===========================================================

function [k,g,s,time] = vdsmex(N,FOV,res,gmax,smax,T,ngmax)

if (nargin < 1)
	error('vdsmex.m:  Need at least 1 argument!');
end;
if (nargin < 2)
	FOV = 24;
end;
if (nargin < 3)
	res = 1;
end;
if (nargin < 4)
	gmax = 3.9;
end;
if (nargin < 5)
	smax = 14500;
end;
if (nargin < 6)
	T = .000004;
end;
if (nargin < 7)
	ngmax = 1000000;	% Don't limit gradient length.
end;

oversamp = 8;		% Factor by which to oversample gradient.
Tg = T/oversamp;	% gradient rate.
Td = T;			% data sampling rate.

krmax = 5/res;		% Max kr.

%	Scale FOV coefficients so that they are normalized.
FOVscale = exp(log(1/krmax)*[0:length(FOV)-1]);
FOV = FOV.*FOVscale;

g = vds_mex(smax,gmax,Tg,Td,N,FOV,krmax,oversamp*ngmax);

g = g(1:oversamp:length(g));
g = [real(g(:)) imag(g(:))];

[k,g,s,m1,m2,time] = calcgradinfo(g,T);



