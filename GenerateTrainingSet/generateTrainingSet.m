% This script generates the training set
% Author: Chiranjib Saha and Harpreet S. Dhillon
% This script generates the sequence of networks and the optimal set of
% simultaneously active links. This script uses the functions in the folder
% ggplab which contains codes for geometric programming [1]

% [1] S. Boyd, S.-J. Kim, L. Vandenberghe, and A. Hassibi, “A tutorial on
% geometric programming,” Optimization and Engineering, vol. 8, no. 1,
%  p. 67, Apr 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;
parameters;
T = 50; % size of the training set
%N = poissrnd(N_av);
for tt = 1:T
 N = poissrnd(N_av);
 [link_distance,tr_loc,rec_loc,S_max,maxrate] =  funGenerateNetwork(N,diskradius);
 %TrainingData.H = H;
 TrainingData.link_distance = link_distance; 
 TrainingData.S_max = S_max; 
 TrainingData.tr_loc = tr_loc;
 TrainingData.rec_loc = rec_loc;
 TrainingData.maxrate = maxrate;
 TrainingDataSet{tt} = TrainingData; 
 fprintf('\ncompleted %f percent',tt/T*100);
end

save('TrainingData','TrainingDataSet','diskRadius'); 