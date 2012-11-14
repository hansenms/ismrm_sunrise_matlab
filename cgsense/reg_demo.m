
load cardiac_data_tutorial;


%Reconstruction with no regularization
figure(1);colormap(gray);
rho_noreg = cg_recon(m,@E_SENSE,e_args);

%Standard Tikhonov regularization
figure(2);colormap(gray);
rho_tikh = cg_recon(m,@E_SENSE,e_args,'fL', @L_std_Tikh,'lambda',0.05,'limit',1e-6);

%Let's combine the coils in the dc_images to use as weights
sum_sq = sum(abs(csm).^2,3);
sum_sq(sum_sq == 0) = max(sum_sq(:)).*1e-10;
dc_img = sum(dc_images .* conj(csm),3) ./ sum_sq;
w = abs(dc_img); clear dc_img;
w(w == 0) = max(w(:)).* 1e-5;
w = w .* (length(w(:))/sum(w(:))); %Normalize weights
e_args.weights = w .^ -1;

figure(3);colormap(gray);
rho_w = cg_recon(m,@E_SENSE,e_args, 'fL', @L_weight,'lambda', 0.1,'limit',1e-6);


%Let's reconstruct a reference image for comparison:
rho_ref = sum(ktoi(data_single_frame,[1,2]) .* conj(csm), 3) ./ sum_sq;

%And show the results
figure(4);colormap(gray);
subplot(2,2,1);
imagesc(abs(rho_ref));axis image; axis off; title('Reference');
cx = caxis;

subplot(2,2,2);
imagesc(abs(rho_noreg));axis image; axis off; title('No regularization');
caxis(cx);

subplot(2,2,3);
imagesc(abs(rho_tikh));axis image; axis off; title('Standard Tikhonov');
caxis(cx);

subplot(2,2,4);
imagesc(abs(rho_w));axis image; axis off; title('Weighted with DC image');
caxis(cx);



