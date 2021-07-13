n_chn_realizations = 20;
Nsub = 64;
Nr = 64;
interference = zeros(Nsub,Nsub,Nr,n_chn_realizations);
sirs = [10,15,20];

dd_in = load('./H_64x64x64_EPA.mat');
H_org = dd_in.noiseless;

for kk = 1:length(sirs)
    sir = sirs(kk);
    fac = sqrt(1/power(10,sir/10));
    disp(fac);
    for i = 1:n_chn_realizations
        x1 = 8*(randperm(8)-ones(1,8));
        x2 = 8*(randperm(8)-ones(1,8));
        H_int = LTE_Eva_Model;
        H_relevant = H_int(1:Nsub,1:8,:);
        H_int2 = repmat(H_relevant, 1, Nsub/8,1);
        x_antenna1 = 2*randi(1,Nsub,Nsub)-1;
        x = repmat(reshape(x_antenna1,[Nsub,Nsub,1]),[1,1,Nr]);
        H_norm_org = reshape(H_org(:,:,:,i),[Nsub*Nr,Nsub]);
        H_norm_int2 = reshape(H_int2,[Nsub*Nr,Nsub]);
        interference(:,:,:,i) = (fac*norm(H_norm_org)/norm(H_norm_int2))*H_int2.*x;
    %     whitening = rand(Nsub)>.95;
        whitening = zeros(Nsub,Nsub);
        for k = 1:2
            x1_k = x1(k);
            x2_k = x2(k);
            whitening(x1_k+1:x1_k+8,x2_k+1:x2_k+8) = ones(8,8);
        end
        whitening_3d = repmat(reshape(whitening,[Nsub,Nsub,1]),[1,1,Nr]);
        interference(:,:,:,i) = interference(:,:,:,i).*whitening_3d;
    end
    save('interference_'+string(sirs(kk))+'.mat','interference');  
end
save('interference_.mat','interference');  