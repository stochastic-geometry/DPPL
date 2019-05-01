% Approximates the MAP of the specified DPP by starting from the empty set
% and iteratively adding the next best element to the selected set.  This
% is the greedy algorithm of Nemhauser and Wolsey (1978).  See readme.txt
% in this folder for a description of inputs and output.
function C = greedyMAP(L, varargin)
nvars = numel(varargin);
assert(nvars == 0 || nvars == 1);

% Greedily add elements one at a time.
C = [];
N = size(L, 1);
U = 1:N;
num_left = N;
while numel(U) > 0
  % Compute the determinant with each remaining unused element added to
  % the chosen set.
  scores = diag(L);
  
  % Select the max-scoring addition to the chosen set.
  [max_score, max_loc] = max(scores);
  if max_score < 1
    break;
  end
  C = [C; U(max_loc)];
  U(max_loc) = [];
  
  % Compute the new kernel, conditioning on the current selection.
  inc_ids = [1:max_loc-1 max_loc+1:num_left];
  L = inv(L + diag([ones(max_loc - 1, 1); 0; ones(num_left - max_loc, 1)]));
  num_left = num_left - 1;
  L = inv(L(inc_ids, inc_ids)) - eye(num_left);
  
  % If enforcing 1:1, throw away any unchosen pairs that would now
  % violate 1:1 if chosen.
  if nvars == 1
    chosen_uv = uv(:, max_loc);
    uv(:, max_loc) = [];
    violators = find((uv(1, :) == chosen_uv(1)) | (uv(2, :) == chosen_uv(2)));
    U(violators) = [];
    uv(:, violators) = [];
    L(violators, :) = [];
    L(:, violators) = [];
    num_left = num_left - numel(violators);
  end
end

% The below code is slow due to the need to compute many determinants;
% while it often improves the return value of this function, it still does
% not allow greedy to beat softmax probability-wise, and comes with a high
% efficiency penalty.
%
% % In the unconstrained case, we can also try greedily deleting elements.
% if nvars == 1
%   return;
% end
% C_up = C;
% C = 1:N;
% prev_score = det(L_saver);
% scores = zeros(N, 1);
% while numel(C) > 0
%   for i = 1:numel(C)
%     proposal = C([1:i-1 i+1:end]);
%     scores(i) = det(L_saver(proposal, proposal));
%   end
%  
%   [max_score, max_loc] = max(scores);
%   if max_score <= prev_score
%     break;
%   end
%   prev_score = max_score;
%   C = C([1:max_loc-1 max_loc+1:end]);
% end
%
% if det(L_saver(C_up, C_up)) > prev_score
%   C = C_up;
% end
