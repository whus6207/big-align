function [ costFrob cost ] = computeCost(A, B, P, Q, lamda, mu)
%Find how good is the alignment between the permuted 
% matrix (A), and matrix B, given the row and column
% permutation matrices P and Q.

costFrob = sum(sum((sparse(P)*A*sparse(Q)-B).^2));
cost = sum(sum((sparse(P)*A*sparse(Q)-B).^2)) + lamda*sum(sum(P)) + mu*sum(sum(Q));


%%%%%%%%%%% CREATE AUGMENTED COST FUNCTION!

end

