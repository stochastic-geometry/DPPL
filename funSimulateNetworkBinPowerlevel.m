function [H,link_distance,tr_loc,rec_loc,S_max,maxrate] = funSimulateNetworkBinPowerlevel(N,diskRadius,N0,alpha,link_dist,mode)
 GenerateNewTopology = 1;
 % % Drop few transreceivers
 %N = 12; % number of links in network
 %diskRadius = 50;
 if GenerateNewTopology ==1
  randangle = 2*pi*rand(N,1);
  randradius  = (diskRadius/7)*sqrt(rand(N,1));
  tr_loc = [randradius.*cos(randangle),randradius.*sin(randangle)];
  randorient = 2*pi*rand(N,1);
  rec_loc = [tr_loc(:,1)+link_dist.*cos(randorient),tr_loc(:,2)+link_dist.*sin(randorient)];
%  TWO clusters
%   tr_loc(1:N/2,1)= tr_loc(1:N/2,1)+diskRadius/2;
%   tr_loc(N/2+1:N,1)= tr_loc(N/2+1:N,1)-diskRadius/2;
%   rec_loc(1:N/2,1)= rec_loc(1:N/2,1)+diskRadius/2;
%   rec_loc(N/2+1:N,1)= rec_loc(N/2+1:N,1)-diskRadius/2;
% Three clusters
%   randangle = 2*pi*rand(N,1);
%   randradius  = (diskRadius/7)*sqrt(rand(N,1));
%   tr_loc = [randradius.*cos(randangle),randradius.*sin(randangle)];
%   randorient = 2*pi*rand(N,1);
%   rec_loc = [tr_loc(:,1)+link_dist.*cos(randorient),tr_loc(:,2)+link_dist.*sin(randorient)];
%   %
%   tr_loc(1:N/3,1)= tr_loc(1:N/3,1)+diskRadius/2;
%   tr_loc(N/3+1:2*N/3,1)= tr_loc(N/3+1:2*N/3,1)-diskRadius/2;
%   tr_loc(2*N/3+1:N,2)= tr_loc(2*N/3+1:N,2)+diskRadius/2;
%   %
%   rec_loc(1:N/3,1)= rec_loc(1:N/3,1)+diskRadius/2;
%   rec_loc(N/3+1:2*N/3,1)= rec_loc(N/3+1:2*N/3,1)-diskRadius/2;
%   rec_loc(2*N/3+1:N,2)= rec_loc(2*N/3+1:N,2)+diskRadius/2;

%  save Topology;
 else
    load Topology;
 end
%  figure;
%  hold on;
%  for linkcount = 1:N
%  patchline([tr_loc(linkcount,1)';rec_loc(linkcount,1)'],[tr_loc(linkcount,2)';rec_loc(linkcount,2)'],'edgecolor','k','edgealpha',1);
%  scatter(tr_loc(linkcount,1),tr_loc(linkcount,2),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','blue');
%  scatter(rec_loc(linkcount,1),rec_loc(linkcount,2),'ro','filled','MarkerFaceAlpha',1,'MarkerFaceColor','red');
%  end
%keyboard;
% figure;
% axis square;
% for linkcount = 1:N
%  patchline([tr_loc(linkcount,1)';rec_loc(linkcount,1)'],[tr_loc(linkcount,2)';rec_loc(linkcount,2)'],'edgecolor','k','edgealpha',1);
%  scatter(tr_loc(linkcount,1),tr_loc(linkcount,2),'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','blue');
%  scatter(rec_loc(linkcount,1),rec_loc(linkcount,2),'ro','filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','red');
% end
% 
% plot([tr_loc(:,1)';rec_loc(:,1)'],[tr_loc(:,2)';rec_loc(:,2)'],'k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel Generation %%
 H = ones(N,N);% randn(N,N)+1j*randn(N,N); % generate complex gaussian channel
 % make H symmetric
 %H  = triu(H);
 %H  = H + transpose(H);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% alpha = 1;
% hold on;
% for linkcount = 1:N
%  patchline([tr_loc(linkcount,1)';rec_loc(linkcount,1)'],[tr_loc(linkcount,2)';rec_loc(linkcount,2)'],'edgecolor','k','edgealpha',0.2);
%  scatter(tr_loc(linkcount,1),tr_loc(linkcount,2),'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','blue');
%  scatter(rec_loc(linkcount,1),rec_loc(linkcount,2),'ro','filled','MarkerFaceAlpha',3/8,'MarkerFaceColor','red');
% end
% box on;
% axis('square');
%H = exprnd(1,N,N);
link_distance= pdist2(rec_loc,tr_loc,'euclidean'); % link_distance(i,j) = ||Tr(j)- Rx(i)|| i.e. Row i corresponds to all distances from Txs to Rec i
if strcmp(mode,'search')
 for iteration = 1: 2^N-1
 S =  decimalToBinaryVector(iteration,N);
 S(S==0) = (0.01); % assign lower transmit power 
%  H_eff = (S'*S).*H;
 % delete zero rows and columns
 %H_eff( all(~abs(H_eff),2), : ) = [];
 %H_eff(:, all(~abs(H_eff),1), : ) = [];
%  link_distance_eff = link_distance;
%  pathloss = link_distance_eff.^(-alpha/2);
%  h = H_eff.*pathloss;
%  Total_power = diag(h*h');
%  serving_power=diag((abs(h).^2));
%  interference = Total_power-serving_power; 
%  %% sanity check
%  if sum(interference<0)&& sum(S==1)>1
%      keyboard;
%  end
%  SINR  = serving_power./(N0+interference);
 P_alloc = S; 
 SINR = funComputesinr(link_distance,P_alloc',N0,alpha);
 %logsumrate(iteration) = sum(log(log2(1+SINR)));
 sumrate(iteration) = sum(log2(1+SINR));
end
% find maximum
%[maxrate,maxindex] = max(logsumrate);
[maxrate,maxindex] = max(sumrate);
 S_max =  decimalToBinaryVector(maxindex,N);
 S_max = S_max>0;
else
 S_max =  [];
 maxrate = [];
end


