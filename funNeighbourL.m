% function L=funNeighbourL(xx,yy,lambda,choiceKernel,sigma,theta,N,M)
% This function file creates an L(-ensemble-)matrix, as detailed in the 
% paper by Blaszczyszyn and Keeler[1](Section IV). 
%
% The quality features or covariates (q_x in the above paper) are based on 
% the nearest neighbour distances. The similiarirty matrix (S in the paper)
% which creates replusion among the points, can be formed from either 
% Gaussian or Cauchy kernel function. 
%
% INPUTS:
% xx and yy are the x and y values of the underlying discrete state space,
% which is usually a realization of a point process on a bounded continuous 
% 2-D space such as a square or a disk.
%
% lambda is the point intensity/density (ie average number of points per
% unit area) of the point process, which is used to rescale the distances.
%
% choiceKernel is a variable that takes value 1 (for Gaussian) or 2 (for
% Cauchy) to select the kernel function. 
%
% sigma is a parameter of kernel function.
%
% theta is a fitting parameters for the quality features/covariates.
%
% N is the number of neighbouring points. 
%
% M is the number of distances between neighbour points. M is optional, but
% when used, M must be equal to zero or N-1.
%
% OUTPUTS: 
% An L-ensemble kernel matrix for a determinantal point process on a 
% discrete space; see [1] for details.
% 
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018.
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.

function [L,SMatrix]=funNeighbourL(Total_power,link_distance,tr_loc,rec_loc,diskradius,choiceKernel,sigma,theta)
 theta=theta(:);

 alpha=1; %an additional parameter for the Cauchy kernel

  %START - Creation of L matrix
 sizeL=size(Total_power,1); %width/height of L (ie cardinality of state space)
 %START -- Create q (ie quality feature/covariage) vector
 % theta = [theta0,theta1,theta2] ?
 %zeroth term
 thetaFeature=theta(1)*ones(sizeL,1);
     
    
    
%     Serving_power =  diag(Total_power); %%%% First feature
%     thetaFeature=thetaFeature ...
%         +sum(theta(2).*Serving_power,2)./10;  %%%%% ADDING some regularization here to make values stable
%     
    dummy = Total_power; 
    dummy(logical(eye(size(Total_power))))=-999;
    dummy_sort =  sort(dummy,2,'descend');
    Best_interference_power =  dummy_sort(1); %%%% check this line % I think the max will operate on rows
    thetaFeature=thetaFeature ...
        +sum(theta(2).*Best_interference_power,2); %%ADDING some regularization here to make values stable
  
   Best_2_interference_power =  dummy_sort(2); %%%% check this line % I think the max will operate on rows
    thetaFeature=thetaFeature ...
        +sum(theta(3).*Best_2_interference_power,2);

qVector=exp(thetaFeature/100); %find q vector (ie feature/covariate values)
%END -- Create q vector

%START - Create similarity matrix S
% if sigma~=0
%     %all squared distances of x/y difference pairs
   %rrDiffSquared=link_distance.^2;
   %rrDiffSquared = rrDiffSquared+ rrDiffSquared';
   xxDiff=bsxfun(@minus,tr_loc(:,1),tr_loc(:,1)'); yyDiff=bsxfun(@minus,tr_loc(:,2),tr_loc(:,2)');
   ttDiffSquared=(xxDiff.^2+yyDiff.^2)./diskradius^2;
   
   xxDiff=bsxfun(@minus,rec_loc(:,1),rec_loc(:,1)'); yyDiff=bsxfun(@minus,rec_loc(:,2),rec_loc(:,2)');
   rrDiffSquared=(xxDiff.^2+yyDiff.^2)./diskradius^2;
    if choiceKernel==1
%         %%Gaussian kernel
         SMatrix=(1./1)*exp(-(ttDiffSquared+rrDiffSquared)/sigma^2);
    elseif choiceKernel==2
%         
%         %%Cauchy kernel
         SMatrix=1./(1+(ttDiffSquared+rrDiffSquared)/sigma^2).^(alpha+1/2);
%     end
    else
     SMatrix=eye(sizeL);    
    end

%   v_mat = Total_power';
%   v_mat_mag = sqrt(sum(v_mat.^2));
%   M = repmat(v_mat_mag,size(v_mat,1),1);
%   v_mat_norm = v_mat./M;
%   SMatrix = v_mat_norm'*v_mat_norm;
%repmat(sqrt(diag(v_mat*v_mat')),1,size(v_mat,2))
% v_mat = v_mat./repmat(sqrt(diag(v_mat*v_mat')),1,size(v_mat,2));
% SMatrix = abs(v_mat*v_mat');  % getting rid of the sign


%SMatrix = 0.5 * (SMatrix + transpose(SMatrix));
%END - Create similarity matrix S with Gaussian kernel

%START Create L matrix
qMatrix=repmat(qVector',size(SMatrix,1),1); %q diagonal matrix
L=(qMatrix').*SMatrix.*qMatrix;
%END Create L matrix

end