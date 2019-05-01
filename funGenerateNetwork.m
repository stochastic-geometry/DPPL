function [link_distance,tr_loc,rec_loc,S_max,maxrate] = funGenerateNetwork(N,diskRadius)
 addpath(genpath('.\ggplab\'));
 tolerance = 0.05;
 beta = 1.1;
 N0= 0.0005;
 alpha = 2;
 link_dist =1;
 % Generate a network topology %
 [~,link_distance,tr_loc,rec_loc,~,~] = funSimulateNetworkBinPowerlevel(N,diskRadius,N0,alpha,link_dist,'no_search');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 P_alloc_guess = ones(N,1);
 SINR_guess = funComputesinr(link_distance,P_alloc_guess,N0,alpha);
 G = link_distance.^-alpha; % this is the pathloss matrix
 G_empty_diag = G;
 G_empty_diag(logical(eye(size(G))))=0;
% encode the variables vector
 loopcount = 0;
 min_var_arr =[];
 while true
  fprintf('\n Solving GP, loop = %d',loopcount);
  loopcount = loopcount+1;
  % note that GP is free from initial guess hyperparameters. 
  gpvar variables(2*N)
  objfunc  = prod(variables(N+1:2*N).^(-(SINR_guess./(1+SINR_guess))));
  constr = [ beta^-1*SINR_guess <= variables(N+1:2*N) ; 
           variables(N+1:2*N) <= beta*SINR_guess ;
           (N0*diag(G).^-1.*variables(1:N).^-1 +...
            (G_empty_diag*variables(1:N)).*diag(G).^-1.*variables(1:N).^-1).*variables(N+1:2*N)<=ones(N,1);
           variables(1:N)<= ones(N,1)];
   
 % solve the power control problem
  [min_var, solution, status] = gpsolve(objfunc, constr);
  assign(solution);       
  SINR_guess = solution{2}(N+1:2*N);    
  min_var_arr(loopcount) = min_var;
  fprintf('\n At iteration %d, min value computed = %f',loopcount,min_var);
  if loopcount>25 && min_var_arr(end-1)-min_var_arr(end)<tolerance
     break
  end
 end
% Sanity_Check %
 P_alloc = solution{2}(1:N);
 %SINR = funComputesinr(link_distance,P_alloc,N0,alpha);  
 S_max = (P_alloc>0.1)';
 P_alloc_quant = (P_alloc<0.1)*0.01+~(P_alloc<0.1)*1;
 SINR = funComputesinr(link_distance,P_alloc,N0,alpha);  
 SINR_quant = funComputesinr(link_distance,P_alloc_quant,N0,alpha); 
 fprintf('\n Sum rate: GP = %f, quant = %f',sum(log2(1+SINR)),sum(log2(1+SINR_quant)));
% 
 maxrate = sum(log2(1+SINR_quant));
% plotNetwork(tr_loc,rec_loc,N,S_max);
end

