clc;
rng('default')
Nsub = 64;
Nr = 64;
SNR_vec = 0:5:20;
n_chn_realizations = 20;
Rhh = zeros(Nsub*Nr);
H_org = zeros(Nsub,Nsub,Nr,n_chn_realizations);
Htotal = zeros(1200,12,Nr,n_chn_realizations);
H_org_noisy = zeros(Nsub,Nsub,Nr,length(SNR_vec),n_chn_realizations);
mse_ls = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse_true = zeros(n_chn_realizations,length(SNR_vec));

for ii = 1:n_chn_realizations
    Htotal(:,:,:,ii) = LTE_Eva_Model;
    H(:,:,:) = Htotal(:,:,:,ii);
    H_relevant = H(1:Nsub,1:8,:);
    H_org(:,:,:,ii) = repmat(H_relevant, 1, Nsub/8,1);
    H_2D = reshape(H_org(1:Nsub,1,:,ii),[Nsub*Nr,1]);
    Rhh = Rhh + (1/n_chn_realizations)*(H_2D*H_2D');
end

for iter = 1:n_chn_realizations
    chn_norm_factor = norm(reshape(H_org(:,1,:,iter),[Nsub,Nr]))/sqrt(Nsub*Nr); 
    noise_matrix = randn(Nsub,Nsub,Nr)+1j*randn(Nsub,Nsub,Nr);

    for k = 1:length(SNR_vec)
        snr = SNR_vec(k);
        std_dev = sqrt(1/power(10,snr/10));
        
        noise = chn_norm_factor * std_dev * noise_matrix; 
        H_org_noisy(:,:,:,k,iter) = H_org(:,:,:,iter) + noise;

        % Pilots are considered as 1 without loss of any generality
        % LS estimation
        H_ls = H_org_noisy(:,1,:,k,iter);
        diff_ls = H_org(1:Nsub,1,:,iter) - H_ls(1:Nsub,1,:);
        mse_ls(iter,k) = norm(reshape(diff_ls,[Nsub,Nr]))^2/norm(reshape(H_org(1:Nsub,1,:,iter),[Nsub,Nr]))^2;
        
        H_LMMSE = Rhh*((Rhh+chn_norm_factor/power(10,snr/10)*eye(Nsub*Nr))\reshape(H_ls,[Nsub*Nr,1]));
        diff_mmse = reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]) - H_LMMSE;
        mse_mmse(iter,k) = norm(diff_mmse)^2/norm(reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]))^2;
        
        H_2D = reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]);
        Rhh_true = H_2D*H_2D';
        H_LMMSE_true = Rhh_true*((Rhh_true+chn_norm_factor/power(10,snr/10)*eye(Nsub*Nr))\reshape(H_ls,[Nsub*Nr,1]));
        diff_mmse_true = reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]) - H_LMMSE_true;
        mse_mmse_true(iter,k) = norm(diff_mmse_true)^2/norm(reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]))^2;
    end

noiseless = H_org;
noisy_0dB = H_org_noisy(:,:,:,1,:);
noisy_5dB = H_org_noisy(:,:,:,2,:);
noisy_10dB = H_org_noisy(:,:,:,3,:);  
noisy_15dB = H_org_noisy(:,:,:,4,:);
noisy_20dB = H_org_noisy(:,:,:,5,:);
end

save('H_64x64x64_EPA.mat','noiseless','noisy_0dB','noisy_5dB','noisy_10dB','noisy_15dB','noisy_20dB');

figure
semilogy(SNR_vec, (mean(mse_ls)), '-+','LineWidth',2)
grid
hold
semilogy(SNR_vec, (mean(mse_mmse)), '-+','LineWidth',2)
semilogy(SNR_vec, (mean(mse_mmse_true)), '-+','LineWidth',2)