n_chn_realizations = 20;
Nsub = 64;  
interference = zeros(Nsub,Nsub,n_chn_realizations);
whitening = zeros(Nsub,Nsub);

dd_in = load('./H_64x64_allSNR.mat');
H_org = dd_in.noiseless;

for i = 1:n_chn_realizations
    H_int = LTE_Eva_Model;    
    H_relevant = H_int(1:Nsub,1:8);
    H_int2 = repmat(H_relevant, 1, Nsub/8);
    x = 2*randi(1,Nsub)-1;
    interference(:,:,i) = (0.5*norm(H_org(:,:,i))/norm(H_int2))*H_int2.*x;
    whitening = rand(Nsub,1)>.90; 
    interference(:,:,i) = interference(:,:,i).*whitening;
end
save('interference.mat','interference');