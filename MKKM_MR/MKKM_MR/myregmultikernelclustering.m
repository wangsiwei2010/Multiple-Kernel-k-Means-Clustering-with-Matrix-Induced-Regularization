function [H_normalized,gamma,obj] = myregmultikernelclustering(K,M,cluster_count,lambda)

nbkernel = size(K,3);
gamma0 = ones(nbkernel,1)/nbkernel;
flag = 1;
iter = 0;
while flag
    iter = iter + 1;
    KC  = mycombFun(K,gamma0);
    H = mykernelkmeans(KC,cluster_count);
    obj(iter)  = calObj(H,K,M,gamma0,lambda);
    [gamma]= updatekernelweightsV2(H,K,M,lambda);
    gamma0 = gamma;
    if iter>2 && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-6)
        flag =0;
    end
end
H_normalized = H./ repmat(sqrt(sum(H.^2, 2)), 1,cluster_count);