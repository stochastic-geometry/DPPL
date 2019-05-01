clear all;close all;
parameters;
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