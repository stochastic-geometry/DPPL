% Author: Chiranjib Saha and Harpreet S. Dhillon
% Main code for "Machine Learning meets Stochastic Geometry:
% Determinantal Subset Selection for Wireless Networks"


% Train the DPP model 
%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
 addpath(genpath('.\DataSet\'));
parameters;
% load the training data
load TrainingDataSingleClusters;

% Set the size of the training set
T =50;
TrainingCollection = {TrainingDataSet{1:T}};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
booleOptSigma = 1; % this option enables the parameterization of S
theta = [1,1,1];
sigma = 1;
param = [theta, sigma]; 
funMax_theta=@(param)funLikelihood_data(T,TrainingCollection,diskradius,choiceKernel,param,alpha);
funMin=@(theta)(-1*funMax_theta(theta)); %define function to be minimized

% Initial values of the parameters
if choiceKernel==3
    thetaGuess = [10,1,1];
else
    thetaGuess = [10,1,1,1];
end

options=optimset('Display','iter'); %options for fminsearch
thetaMax=fminunc(funMin,thetaGuess,options); %minimize function 

if booleOptSigma
    sigma=thetaMax(end); %retrive sigma values from theta
    thetaMax=thetaMax(1:end-1);
end
% Display the trained parameters
sigma
thetaMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of DPP training%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% learn the activation probability
p_a_1= 0; 
p_a_2 = 0;
for tt = 1: T
 p_a_1 = p_a_1+ TrainingCollection{tt}.S_max(1);
 p_a_2 =  p_a_2+TrainingCollection{tt}.S_max(end);
end
p_a_1 = p_a_1/T;
p_a_2 = p_a_2/T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Start of Testing %%
load('TrainingDataSingleClustersPart2.mat');
TestSet =  TrainingDataSet;
numbSim =size(TestSet,2);
for ss=1:numbSim
    % Read the dataset %%
    link_distance = TestSet{ss}.link_distance; 
    tr_loc = TestSet{ss}.tr_loc; 
    rec_loc = TestSet{ss}.rec_loc; 
    S_max = TestSet{ss}.S_max; 
    maxrate = TestSet{ss}.maxrate;
    N = size(TestSet{ss}.link_distance,2);
    %%%%%%%%%%%%%%%%%%%%%%
    Total_power= (link_distance.^(-alpha/2)).^2;
    [L,~] = funNeighbourL(Total_power,link_distance,tr_loc,rec_loc,diskradius,choiceKernel,sigma,thetaMax);
    S = sample_dpp(decompose_kernel(L));
    
    P_alloc_DPP_sample = 0.01*ones(1,N);
    P_alloc_DPP_sample(S) = 1;
    

    SINR = funComputesinr(link_distance,P_alloc_DPP_sample',N0,alpha);
    sumrate_Random_DPP(ss) = sum(log2(1+SINR)); 
    
    S_MAP = greedy_sym(L);
    P_alloc_MAP = 0.01*ones(1,N);
    P_alloc_MAP(S_MAP) = 1;
    SINR = funComputesinr(link_distance,P_alloc_MAP',N0,alpha);
    sumrate_MAP(ss) = sum(log2(1+SINR));
    
    S = rand(1,N)<p_a_1;
    P_alloc_random = (S==0)*0.01+(S==1)*1;
    SINR = funComputesinr(link_distance,P_alloc_random',N0,alpha);
    sumrate_random(ss) = sum(log2(1+SINR)); 
    sumrate_opt(ss) = maxrate;
    
    fprintf('\n Sum rate: optimal = %f, DPP = %f, Random = %f',sumrate_opt(ss),sumrate_MAP(ss),sumrate_random(ss));
end
save('rate_compare_data','sumrate_opt','sumrate_MAP','sumrate_random','sumrate_Random_DPP')

% plotting module %
close all
figure(1);
axes1 = axes('Parent',figure(1));
hold(axes1,'on');
[f,x] =ecdf(sumrate_opt);
l_opt = plot(x,f,'-','linewidth',2,'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);
[f,x] =ecdf(sumrate_MAP);
l_DPP_MAP = plot(x,f,'k--','linewidth',2);
[f,x] =ecdf(sumrate_random);
l_ind = plot(x,f,'linewidth',2,'LineStyle','-.',...
    'Color',[0.231372549019608 0.443137254901961 0.337254901960784])
[f,x] = ecdf(sumrate_Random_DPP);
l_DPP_sample = plot(x,f,'r-','linewidth',2);
l =legend([l_opt,l_DPP_MAP,l_DPP_sample,l_ind],'Optimum','DPP (MAP inference)','DPP (Sampling)','Independent Thinning');
set(l,'interpreter','latex','fontsize',16);
xlabel('Sum rate (bps)','interpreter','latex','fontsize',16);
ylabel('CDF','interpreter','latex','fontsize',16);
box on;
grid on;
xlim([8,27])
set(axes1,'FontName','Times New Roman','FontSize',14);