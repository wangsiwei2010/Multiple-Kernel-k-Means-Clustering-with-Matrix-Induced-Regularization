function [gamma,obj]= updatekernelweightsV2(T,K,HE0,lambda)

num = size(K,1);
nbkernel = size(K,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0 = eye(num)-T*T';
ZH = zeros(nbkernel,1);
for p = 1 : nbkernel
    ZH(p) = trace(K(:,:,p)*U0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = lambda*HE0+2*diag(ZH);
f = zeros(nbkernel,1);
A = [];
b = [];
Aeq = ones(nbkernel,1)';
beq = 1;
lb  = zeros(nbkernel,1);
ub =  ones(nbkernel,1);

[gamma,obj]= quadprog(H,f,A,b,Aeq,beq,lb,ub);
gamma(gamma<1e-6)=0;
gamma = gamma/sum(gamma);