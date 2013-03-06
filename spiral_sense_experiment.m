%clear
clear all;

FOV = 30; %cm
SLEW = 14000; %gauss per cm per s
GRAD = 2.4; %G/cm
Ts = 5e-6; %Sample time s
Nint = 12;
MATRIX = 256;
KMAX = (1/(FOV/MATRIX))/2;
GRADMAX = 1e4;
Acc = 4;

trajectory = 'cartesian';

if (~strcmp(trajectory,'cartesian')),
    %Calculate trajectory and weights
    if (1),
        %Spiral trajectory
        [g,k,w] = spiral_trajectory(SLEW,GRAD,Ts,Ts,Nint,FOV,KMAX,GRADMAX);
        k = (k ./ (2*KMAX))*MATRIX;
    else
        %Radial trajectory
        Nint = MATRIX;
        kx = -bitshift(MATRIX,-1):1:(bitshift(MATRIX,-1)-1);
        k = zeros(length(kx),Nint);
        for ii=1:Nint,
           an = ((ii-1)/(Nint-1))*pi;
           k(:,ii) = complex(cos(an)*kx,sin(an)*kx);
           w(:,ii) = abs(kx);
        end
        w(w == 0) = 0.25;
        w = w(:);
        k = k(:);
    end
end


%Load Data
load ../im1.mat
load ../smaps_phantom.mat

%Coil data in k-space on uniform grid
coil_data = ismrm_transform_image_to_kspace(repmat(im1,[1 1 size(smaps,3)]) .* smaps,[1 2]);

%Sample on trajectory
if (~strcmp(trajectory,'cartesian')),
    for c=1:size(smaps,3), spiral_coil_data(:,c) = grid_data_bck([real(k) imag(k)], coil_data(:,:,c),[0 0]); end;
    noise_level = max(abs(im1(:)))*0.05;
else
    noise_level = 0.01*sqrt(sum(abs(coil_data(:)).^2)/numel(coil_data));
end

%Adding some noise
%noise_level = 0.001*sqrt(sum(abs(spiral_coil_data(:)).^2)/numel(spiral_coil_data));
%noise = complex(randn(size(spiral_coil_data)),randn(size(spiral_coil_data))).*noise_level; 
%spiral_coil_data = spiral_coil_data + noise;


if (strcmp(trajectory,'cartesian')),
    noise = complex(randn(size(coil_data)),randn(size(coil_data))).*noise_level;
    grid_coil_data_sub = zeros([size(coil_data) Acc]);
    for a=1:Acc,
        for c=1:size(smaps,3), grid_coil_data_sub(:,a:Acc:end,c,a) = (coil_data(:,a:Acc:end,c) + noise(:,a:Acc:end,c)) * Acc; end    
        sub_coil_images(:,:,:,a) = ismrm_transform_kspace_to_image(grid_coil_data_sub(:,:,:,a),[1,2]);
    end
else
    noise = complex(randn(size(spiral_coil_data)),randn(size(spiral_coil_data))).*noise_level;

    %Accelerated date
    noise = reshape(noise,size(noise,1)/Nint,Nint,size(smaps,3));

    spiral_coil_data = reshape(spiral_coil_data,size(spiral_coil_data,1)/Nint,Nint,size(smaps,3));
    w = reshape(w,size(w,1)/Nint,Nint);
    k = reshape(k,size(k,1)/Nint,Nint);

    for a=1:Acc,
        for c=1:size(smaps,3), grid_coil_data_sub(:,:,c,a) = grid_data([real(vec(k(:,a:Acc:end))) vec(imag(k(:,a:Acc:end)))],spiral_coil_data(:,a:Acc:end,c)+noise(:,a:Acc:end,c),vec(w(:,a:Acc:end)),[MATRIX MATRIX],[0 0]); end;
        sub_coil_images(:,:,:,a) = ismrm_transform_kspace_to_image(grid_coil_data_sub(:,:,:,a),[1,2]);
    end
end

%Reconstruct fully sampled data
%for c=1:size(smaps,3), grid_coil_data(:,:,c) = grid_data([real(k) imag(k)], spiral_coil_data(:,c)+noise(:,c),w(:),[MATRIX MATRIX],[0 0]); end;

full_coil_images = sum(sub_coil_images,4)/(Acc*Acc);%ismrm_transform_kspace_to_image(grid_coil_data,[1,2]);
csm = ismrm_estimate_csm_inati(full_coil_images);
acomb_image = sum(conj(csm).*full_coil_images,3);





%Calculate SENSE weights
% adaptively combined weights with aliasing control
Nx = MATRIX;
Ny = MATRIX;
Nc = size(smaps,3);

Wb = zeros(Nx,Ny,Nc);
for x = 1:Nx
    for y = 1:Ny      
        d = permute(sub_coil_images(x,y,:,:), [4 3 2 1]);
        size(d)
        t = acomb_image(x,y) * ones(Acc,1);
        size(d)
        size(t)
        v = d\t;
        size(v)
        break
        Wb(x,y,:) = v;
    end
    break
end


iters = 100;
for it = 1:iters,
fprintf('Iteration %d/%d\n',it, iters);

%Generate "new" sub_coil images (with different noise)
if (strcmp(trajectory,'cartesian')),
    noise = complex(randn(size(coil_data)),randn(size(coil_data))).*noise_level;
    grid_coil_data_sub = zeros([size(coil_data) Acc]);
    for a=1:Acc,
        for c=1:size(smaps,3), grid_coil_data_sub(:,a:Acc:end,c,a) = (coil_data(:,a:Acc:end,c) + noise(:,a:Acc:end,c)) * Acc; end    
        sub_coil_images(:,:,:,a) = ismrm_transform_kspace_to_image(grid_coil_data_sub(:,:,:,a),[1,2]);
    end
else
    noise = complex(randn(size(spiral_coil_data)),randn(size(spiral_coil_data))).*noise_level;
    for a=1:Acc,
        for c=1:size(smaps,3), grid_coil_data_sub(:,:,c,a) = grid_data([real(vec(k(:,a:Acc:end))) vec(imag(k(:,a:Acc:end)))],spiral_coil_data(:,a:Acc:end,c)+noise(:,a:Acc:end,c),vec(w(:,a:Acc:end)),[MATRIX MATRIX],[0 0]); end;
        sub_coil_images(:,:,:,a) = ismrm_transform_kspace_to_image(grid_coil_data_sub(:,:,:,a),[1,2]);
    end
end
%Reconstruct image:
for a=1:Acc,
    I(:,:,a,it) = sum(Wb.*sub_coil_images(:,:,:,a),3);
end

end
