% Approximates the MAP of the specified DPP by starting from the empty set
% and iterating *once* through the items, adding as specified in Buchbinder
% et al. (2012).  See readme.txt in this folder for a description of inputs
% and output.
function S = greedy_sym(L)
% Initialize.
N = size(L,1);
S = [];

% Iterate through items.
for i = 1:N
  lift0 = obj(L,[S i]) - obj(L,S);
  lift1 = obj(L,[S (i+1):N]) - obj(L,[S i:N]);
  
  if lift0 > lift1
    S = [S i];
  end
end

function f = obj(L, S)
f = log(det(L(S, S)));