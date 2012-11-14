function [unmix, b1, unmix_sc] = ismrm_calculate_grappa_unmixing(full_data, kernel_size, acc_factor, target_data, data_mask)
%
%   INPUT:
%       full_data [read,phase,coil]                  : Fully sampled k-space
%       kernel_size [freq_positions, phase_positions]: e.g. [4 5]
%       acc_factor                                   : Acceleration factor
%
%   OUTPUT:
%       unmix [x,y,coil]                             : Image unmixing coefficients
%
%   Typical usage:
%       [unmix] = calculate_grappa_unmixing(time_average, [5 4], 4);
%
%

if (nargin<4),
    target_data = full_data;
else
    if (isempty(target_data)),
        target_data = full_data;
    end
end

if (nargin<5),
    data_mask = ones(size(full_data,1),size(full_data,2));
end

if (length(size(full_data)) == 2),
    coils = 1;
else
    coils = size(full_data,length(size(full_data)));
end

if (length(size(target_data)) == 2),
    target_coils = 1;
else
    target_coils = size(target_data,length(size(target_data)));
end

coefficients = kernel_size(1)*kernel_size(2)*coils
overdetermined_factor = 10;

[d1_min,d2_min] = ind2sub(size(data_mask),find(data_mask,1,'first'));
[d1_max,d2_max] = ind2sub(size(data_mask),find(data_mask,1,'last'));

d1_range = (bitshift(kernel_size(1),-1)+d1_min):(d1_max-bitshift(kernel_size(1)+1,-1));
d2_range = (bitshift(kernel_size(2)*acc_factor,-1)+d2_min):(d2_max-bitshift(kernel_size(2)*acc_factor+1,-1));

k_locations = length(d1_range)*length(d2_range)

%Let's reduce the number of equations if we are too overdetermined
if (0),
if (k_locations > coefficients*overdetermined_factor),
    if (length(d2_range) > sqrt(coefficients*overdetermined_factor));
        d2_locations = ceil(sqrt(coefficients*overdetermined_factor))
        d1_locations = ceil(sqrt(coefficients*overdetermined_factor))
    else,
        d2_locations = length(d2_range)
        d1_locations = ceil((coefficients*overdetermined_factor)/d2_locations)
    end
    d1_range = d1_range([1:d1_locations]+bitshift(length(d1_range)-d1_locations,-1))
    d2_range = d2_range([1:d2_locations]+bitshift(length(d2_range)-d2_locations,-1))
    k_locations = length(d1_range)*length(d2_range)
end
end

kernel = zeros(kernel_size(1),kernel_size(2)*acc_factor,coils,target_coils);
fprintf('Calculating grappa kernels...\n');
for s=1:(acc_factor),
    fprintf('Inversions %d of %d...', s, (acc_factor));
    A = zeros(k_locations,coefficients);
    b = zeros(k_locations,target_coils);
    
    k_loc_counter = 1;
    for d1=d1_range,
        d1_vals = [d1:d1+kernel_size(1)-1]-bitshift(kernel_size(1),-1);
        for d2=d2_range,
            %d2_vals = d2+(([0:(kernel_size(2)-1)]*acc_factor)+(acc_factor-s+1)-bitshift(size(kernel,2),-1)-1);
            d2_vals = d2+(([0:(kernel_size(2)-1)]*acc_factor)+(s)-bitshift(size(kernel,2),-1)-1)+1;
%             k_loc_counter
%             d1_vals
%             d2_vals
%             A(k_loc_counter,:)
%             full_data(d1_vals,d2_vals,:)
            A(k_loc_counter,:) = vec(full_data(d1_vals,d2_vals,:));
            b(k_loc_counter,:) = target_data(d1,d2,:);
            k_loc_counter = k_loc_counter + 1;
        end
    end
    
    fprintf('inverting...');
    if (0),
        %No regularization
        A_inv = pinv( A'*A)*A';
    elseif(0),
        %Truncated SVD
        [U,S,V] = svd(A'*A,0);
        S_inv = zeros(size(S));
        idx = find(S > (1e-4*max(S(:))));
        S_inv(idx) = S(idx).^-1;
        S_inv = S_inv';
        A_inv = (V*S_inv*U')*A';
    else,
        %Tikhonov
        S = svd(A,0);
        A_inv = pinv(A'*A + eye(size(A'*A)).*(1e-3*max(abs(S(:)))).^2)*A';
    end
    
    kernel_set = A_inv*b;
    for c=1:target_coils,
        %kernel_set = A_inv*b(:,c);
        %kernel(:,([0:(kernel_size(2)-1)]*acc_factor)+(acc_factor-s+1),:,c) = reshape(kernel_set,kernel_size(1),kernel_size(2),coils);
        kernel(:,([0:(kernel_size(2)-1)]*acc_factor)+(s+1),:,c) = reshape(kernel_set(:,c),kernel_size(1),kernel_size(2),coils);
    end
    fprintf('done.\n')
end


%kernel(:) = 0;
%for c=1:coils,
%    kernel(bitshift(size(kernel,1),-1)+1,bitshift(size(kernel,2),-1)+1, c, c) = 1;
%end

kernel = flipdim(flipdim(kernel,1),2); %Flip dimensions in preparation for convolution.


fprintf('Calculating B1 map...');
filter_x = hamming(size(data_mask,1));
filter_y = zeros(1,size(data_mask,2));
filter_y(data_mask(1,:)>0) = hamming(sum(data_mask(1,:)));
raw_filter = filter_x*filter_y;

if (target_coils == 1),
    b1 = ones(size(target_data));
else
    b1 = calculateB1map(ktoi(target_data.*repmat(raw_filter,[1 1 target_coils]),[1,2])/sqrt(size(target_data,1)*size(target_data,2)));
end
fprintf('done.\n');

unmix = zeros(size(full_data));
if (nargout > 2),
   unmix_sc = zeros(size(unmix,1),size(unmix,2),coils,coils); 
end
fprintf('Doing B1 weighted combination....');
for c=1:target_coils,
    kernel_pad = ktoi(pad_grappa_kernel(kernel(:,:,:,c),size(target_data)),[1,2]);
    kernel_pad = kernel_pad*(size(kernel_pad,1)*size(kernel_pad,2)/acc_factor);
    if (nargout > 2),
        unmix_sc(:,:,:,c) = kernel_pad;
    end
    unmix = unmix + (kernel_pad .* repmat(conj(b1(:,:,c)),[1 1 coils])); %Should divide by SS of coils, but it is '1' by definition here
end
unmix = unmix;

fprintf('done.\n');

return


function padded_kernel = pad_grappa_kernel(gkernel, image_size)
    padded_kernel = zeros(image_size(1),image_size(2),size(gkernel,3));
    padded_kernel([1:size(gkernel,1)]+bitshift(image_size(1)-size(gkernel,1)-1,-1)+1, ...
        [1:size(gkernel,2)]+bitshift(image_size(2)-size(gkernel,2)-1,-1)+1, :) = gkernel;
return
        
function [Data] = ktoi(Data, Dims)
    if nargin<2,
        Data = fftshift(ifftn(ifftshift(Data)));
    else
        for d=[1:length(Dims)],
            Data = fftshift(ifft(ifftshift(Data),[],Dims(d)));
        end
    end
return
