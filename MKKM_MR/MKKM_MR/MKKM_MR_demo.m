clear
clc
warning off;

path = '.\';
addpath(genpath(path));
dataName = 'proteinFold'; %%% flower17; flower102; proteinFold,caltech101_mit,UCI_DIGIT,ccv
%% %% washington; wisconsin; texas; cornell
%% caltech101_nTrain5_48
%% proteinFold
load([path,'datasets\',dataName,'_Kmatrix.mat'],'KH','Y');
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
qnorm = 2;

% % %%---MKKM-MR(AAAI2016)---%%
tic
Htrntrn = calM(KH);
lambdaset8 = 2.^[-15:1:10];
accval8 = zeros(length(lambdaset8),1);
nmival8 = zeros(length(lambdaset8),1);
purval8 = zeros(length(lambdaset8),1);
for il =1:length(lambdaset8)
    [H_normalized8,gamma81,obj81] = myregmultikernelclustering(KH,Htrntrn,numclass,lambdaset8(il));
    res8 = myNMIACC(H_normalized8,Y,numclass);
    accval8(il) = res8(1);
    nmival8(il) = res8(2);
    purval8(il) = res8(3);
end
t = toc;
