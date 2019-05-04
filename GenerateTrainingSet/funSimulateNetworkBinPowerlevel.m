function [H,link_distance,tr_loc,rec_loc,S_max,maxrate] = funSimulateNetworkBinPowerlevel(N,diskRadius,link_dist)
% Author: Chiranjib Saha and Harpreet S. Dhillon
% This function drops N links on a circular window of diskRadius 
% Output: H : NXN matrix with each entry being 1. This is a dummy output.
% One can implement any fading gain here. 
%         link_distance : NxN matrix with link distances 
%         tr_loc , rec_loc : Nx2 matrix with tx and rx locations
%         S_max : empty matrix
%         maxrate : empty matix

  
  randangle = 2*pi*rand(N,1);
  randradius  = (diskRadius/7)*sqrt(rand(N,1));
  tr_loc = [randradius.*cos(randangle),randradius.*sin(randangle)];
  randorient = 2*pi*rand(N,1);
  rec_loc = [tr_loc(:,1)+link_dist.*cos(randorient),tr_loc(:,2)+link_dist.*sin(randorient)];

%%%%%%%%%%%%%% Code dump: can be uncommented for plotting purposes
%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
S_max =  [];
maxrate = [];


