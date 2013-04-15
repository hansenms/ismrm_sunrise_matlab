%
%	Examples using vds.m and vdsmex.m
%

%	This is a 16-interleave spiral that achieves
%	1mm resolution, and density decreases linearly
%	from supporting 24cm FOV at |k|=0 to 12cm FOV
%	at |k|=maximum.
%

smax = 15000;	 % 150 T/m/s
gmax = 4;	 % G/cm
T = .000004;	 % Seconds
N = 16;		 % Interleaves
Fcoeff = [24 -12] 	% FOV decreases linearly from 24 to 12cm.
res = 1;
rmax = 5/res;		% cm^(-1), corresponds to 1mm resolution.

disp('Calculating Gradient');
[k,g,s,time,r,theta] = vds(smax,gmax,T,N,Fcoeff,rmax);

disp('Plotting Gradient');
g = [real(g(:)) imag(g(:))];
plotgradinfo(g,T);


%
%	Here the example is repeated with vdsmex, which
%	should be a lot faster!

disp('Press Enter to repeat with vdsmex.');
pause;
figure;

ngmax = 1000000;	% Alternative stopping condition,
			% if duration is to be limited.

disp('Calculating Gradient');
[k,g,s,time] = vdsmex(N,Fcoeff,res,gmax,smax,T,ngmax);

disp('Plotting Gradient');
plotgradinfo(g,T);


