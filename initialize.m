function [P, Q, P0, Q0] = initialize(A, B)
% Initial matchings between the nodes with the same degree. If there are
% more than 2 nodes with the same degree, we initialize all the possible 
% matchings uniformly.
% Inputs:
%   A: p2xq1 adjacency matrix of 1st graph
%   B: p1xq2 adjacency matrix of 2nd graph
% Outputs:
%   P: p1xp2 correspondence matrix between nodes in the 1st set of A, B
%   Q: q1xq2 correspondence matrix between nodes in the 2nd set of A, B


% dimensions of correspondence matrices P and Q
p1 = size(B,1);
p2 = size(A,1);
q1 = size(A,2);
q2 = size(B,2);

% degree column-vector of "individual" nodes
D_A1 = full(sum(A,2));
D_B1 = sum(B,2);
% degree row-vector of "community" nodes
D_A2 = full(sum(A,1));
D_B2 = sum(B,1);

% First, initialize the "new" correspondence matrices.
P = zeros(p1, p2);
Q = zeros(q1, q2);

% Initialize P matrix: if node 1 of graph A corresponds to node 3 of graph
% B, then P_31 = 1.
% Trying the following initialization: 
%   Each node, b, in the second graph is matched to each node, a, in the first
%   graph with probability min(degree_a, degree_b)/max(degree_a, degree_b)

Drep_A1 = repmat(D_A1', p1, 1);
Drep_B1 = repmat(D_B1, 1, p2);
degRatioP = Drep_B1./Drep_A1; 
P0 = min(degRatioP, 1./degRatioP);
P0(P0<1)=0;
[idxI idxJ] = find(P0);
for i = 1 : length(idxI)
        P0(idxI(i),idxJ(i)) = P0(idxI(i),idxJ(i))/(nnz(D_A1==D_A1(idxJ(i)))*nnz(D_B1==D_B1(idxI(i))));
end
%nonZerosPerRowP = full(sum((P0~=0),2));
%nonZerosPerRowP_rep = repmat(nonZerosPerRowP, 1, p2); 
%P0 = P0 ./ nonZerosPerRowP_rep;

Drep_A2 = repmat(D_A2', 1, q2);
Drep_B2 = repmat(D_B2, q1, 1);
degRatioQ = Drep_B2./Drep_A2; 
Q0 = min(degRatioQ, 1./degRatioQ);
Q0(Q0<1)=0;
[idxIq idxJq] = find(Q0);
% initialize each matching to 1/(possible node matches with the same degree)
for i = 1 : length(idxIq)
        Q0(idxIq(i),idxJq(i)) = Q0(idxIq(i),idxJq(i))/(nnz(D_A2==D_A2(idxIq(i)))*nnz(D_B2==D_B2(idxJq(i))));
end
%nonZerosPerRowQ = full(sum((Q0~=0),2));
%nonZerosPerRowQ_rep = repmat(nonZerosPerRowQ, 1, q2); 
%Q0 = Q0 ./ nonZerosPerRowQ_rep;



% [vA1 idxA1] = sort(D_A1);
% [vB1 idxB1] = sort(D_B1);
% [vA2 idxA2] = sort(D_A2);
% [vB2 idxB2] = sort(D_B2);
% % Only match nodes of Set 1 that have degree above the average.
% avg_deg1 = min(mean(D_A1), mean(D_B1));
% idxA1_toMatch = idxA1(vA1 > avg_deg1);
% idxB1_toMatch = idxB1(vA1 > avg_deg1);
% % Only match nodes of Set 2 that have degree above the average.
% avg_deg2 = min(mean(D_A2), mean(D_B2));
% idxA2_toMatch = idxA2(vA2 > avg_deg2);
% idxB2_toMatch = idxB2(vA2 > avg_deg2);
% % First, initialize everything to uniform
% P0 = ones(pDim1, pDim2) ./ (pDim2-length(idxA1_toMatch));
% Q0 = ones(qDim1, qDim2) ./ (qDim2-length(idxA2_toMatch));
% % Then, initialize cleverly the top k nodes
% % Find k by thresholding on the degree - match nodes with degree above average.
% % If multiple nodes have the same degree, then initialize the
% % matching for all of them.
% P0(idxB1_toMatch, :) =0;
% P0(:, idxA1_toMatch) =0;
% for i = 1 : length(idxA1_toMatch)
%     P0(idxB1_toMatch(i), idxA1_toMatch(i)) = 0.9;
% end
% Q0(idxA2_toMatch, :) =0;
% Q0(:, idxB2_toMatch) =0;
% for i = 1 : length(idxA2_toMatch)
%     Q0(idxA2_toMatch, idxB2_toMatch) = 0.9;
% end
% % Finally, initialize the new permutation matrices.
% P = zeros(pDim1, pDim2);
% Q = zeros(qDim1, qDim2);

end

