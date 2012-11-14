function [Data] = itok(Data, Dims)
%
% [Data] = itok(Data, Dims)
% 

if nargin<2,
  Data = fftshift(fftn(ifftshift(Data)));
else
  for d=[1:length(Dims)],
    Data = fftshift(fft(ifftshift(Data),[],Dims(d)));
  end
end

return

