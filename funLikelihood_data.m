function logLikelihood=funLikelihood_data(T,TrainingDataSet,diskradius,choiceKernel,param,alpha)
% Code imported from the repository : https://github.com/hpaulkeeler/DetPoisson_MATLAB




 logLikelihoodVector=zeros(T,1);

 

 
%Loop through all training/learning samples
for tt=1:T        
    %Create L matrix (ie for Phi) based on nearest neighbours
    Total_power = ( TrainingDataSet{tt}.link_distance.^(-alpha/2)).^2;
    tr_loc = TrainingDataSet{tt}.tr_loc;
    rec_loc = TrainingDataSet{tt}.rec_loc;
    [L,S]=funNeighbourL(Total_power,TrainingDataSet{tt}.link_distance,tr_loc,rec_loc,diskradius,choiceKernel,param(end),param(1:end-1));
    %Create sub L matrix (ie for Psi)
    index =  find(TrainingDataSet{tt}.S_max>0); 
    subL=L(index,index);
    logLikelihoodVector(tt)=(log(det(subL))-log(det(L+eye(size(L)))));
end
logLikelihood=sum(logLikelihoodVector);
end