function obj = calObj(T,K,H0,gamma0,lambda)

nbkernel = size(K,3);
num = size(K,1);
f = zeros(nbkernel,1);
for p = 1 : nbkernel
    f(p) = (trace(K(:,:,p))-trace(T'*K(:,:,p)*T));
end

obj = (1/2)*gamma0'*(lambda*H0+2*diag(f))*gamma0;
