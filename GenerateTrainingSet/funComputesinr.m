function SINR = funComputesinr(link_distance,P_alloc,N0,alpha)
% Author: Chiranjib Saha and Harpreet S. Dhillon
% This function computes the SINR seen by each receiver in the network. 
% Input: link_distance : NxN matrix
%        P_alloc: power allocation for each link
%        N0 : noise power
%        alpha: pathloss
 link_distance_eff = link_distance;
 pathloss = link_distance_eff.^(-alpha/2);
 Power = repmat(P_alloc',length(P_alloc),1);
 h = sqrt(Power).*pathloss;
 Total_power = diag(h*h');
 serving_power=diag((abs(h).^2));
 interference = Total_power-serving_power; 
 %% sanity check
 if sum(interference<0)&& sum(S==1)>1
     keyboard;
 end
 SINR  = serving_power./(N0+interference);
 %SINR = SINR';
end