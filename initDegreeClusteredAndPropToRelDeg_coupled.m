function [ P, Q, P0, Q0, clustersNoSet1, clustersNoSet2, l_A1, l_A2, l_B1, l_B2 ] = ...
    initDegreeClusteredAndPropToRelDeg_coupled( A, B, ~ )
%% CLUSTERED initialization of matrices P0 and Q0
% Choose how many clusters we want to initialize uniformly
% Inputs:
%   A: p2xq1 adjacency matrix of 1st graph
%   B: p1xq2 adjacency matrix of 2nd graph
%   clustersNo: number of clusters to initialize uniformly
% Outputs:
%   P: p1xp2 correspondence matrix between nodes in the 1st set of A, B
%   Q: q1xq2 correspondence matrix between nodes in the 2nd set of A, B
%
%
% Assumption: Set 2 gives some first matchings between high-degree nodes.
%             Use these matchings in order to initialize P0 (coupled
%             alignment) by using info about the neighbors.
%
%  Algorithm Description:
%  * STEP 1: using an idea similar to scree plot, we find up to which
%       degree we are going to match the top nodes 1-by-1 in the two
%       graphs. We find the knee in the unique degrees of the two graphs,
%       get the minimum rank of the knees and match the nodes up to this
%       point 1-by-1.
%  * STEP 2: we cluster the rest nodes in clustersNo (given by the user)
%       and initialize uniformly these clusters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
%% dimensions of correspondence matrices P and Q
p1 = size(B,1);
p2 = size(A,1);
q1 = size(A,2);
q2 = size(B,2);

%% degree column-vector of "individual" nodes
% for weighted matrices
Bceil = ceil(B);
% for unweighted matrices
D_A1 = full(sum(A,2)');
D_B1 = sum(Bceil,2);
%% degree row-vector of "community" nodes
D_A2 = full(sum(A,1))';
D_B2 = sum(Bceil,1);

%% First, initialize the "new" correspondence matrices.
P = zeros(p1, p2);
Q = zeros(q1, q2);

P0 = zeros(p1, p2);
Q0 = zeros(q1, q2);

%% Nodes of set 1 matched during the alignment of set 2 (coupling)
idx_A1_matched = [];
idx_B1_matched = [];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% sort unique degrees of set 2 in descending order

A2_deg_unique = sort(unique(D_A2), 'descend');
B2_deg_unique = sort(unique(D_B2), 'descend');
l_A2 = length(A2_deg_unique);
l_B2 = length(B2_deg_unique);


%% Prepare the initialization of P0 based on the neighborhoods of set2's nodes
degArep = repmat(D_A1, length(D_A1), 1);
degBrep = repmat(D_B1, 1, length(D_B1));
relDiffP = 2*abs(degArep - degBrep)./(degArep + degBrep);
relDiffP(isnan(relDiffP)) = 0;
%denomP = 1./(1+relDiffP);

%% Prepare the initialization of Q0 so that within the clusters we initialize
% the nodes based on their degrees
degArep = repmat(D_A2, 1,length(D_A2));
degBrep = repmat(D_B2, length(D_B2), 1);
relDiffQ = 2*abs(degArep - degBrep)./(degArep + degBrep);
relDiffQ(isnan(relDiffQ)) = 0;
%denomP = 1./(1+relDiffQ);

%% match the high degree nodes (up to the knee in the degrees plot)
knee1 = detectKneeAtScreeLikePlot( A2_deg_unique );
knee2 = detectKneeAtScreeLikePlot( B2_deg_unique );
threshold2 = min( knee1, knee2 );
for i = 1 : threshold2
    idx_A2 = find(D_A2 == A2_deg_unique(i));
    idx_B2 = find(D_B2 == B2_deg_unique(i));
    total_elements = length(idx_A2) * length(idx_B2);
    Q0(idx_A2, idx_B2) = 1 / total_elements;
    [neighborsA,~] = find(A(:,idx_A2));
    [neighborsB,~] = find(B(:,idx_B2));
    if A2_deg_unique(i) ~= size(A,1) && B2_deg_unique(i) ~= size(B,1)
        neighborsA_updated = neighborsA( ~ismember(neighborsA, idx_A1_matched ) );
        neighborsB_updated = neighborsB( ~ismember(neighborsB, idx_B1_matched ) );
        % Match the neighbors with prob 1/log(difference in their degrees).
        % If (difference in their degrees) = 0, set the probability to 1.
        P0(neighborsA_updated, neighborsB_updated) = 1 ./ ( 1 + relDiffP(neighborsA_updated, neighborsB_updated) );
        idx_A1_matched = [ idx_A1_matched; neighborsA_updated ];
        idx_B1_matched = [ idx_B1_matched; neighborsB_updated ];
    end
end

%% cluster the rest of the degrees in clustersNo and initialize uniformly
%% each cluster
A2_rest = A2_deg_unique(threshold2+1:l_A2);
B2_rest = B2_deg_unique(threshold2+1:l_B2);
% lA2_rest = length(A2_rest);
% lB2_rest = length(B2_rest);
% minlength = min(lA2_rest, lB2_rest);
% clustersNo = min( clustersNo, minlength );
% clustersNoSet2 = clustersNo;
if length(A2_rest) == 1
    clusterMemberships1 = 1;
else
    clusterMemberships1 = clusterdata(A2_rest, 1);
end
clustersNo = length(unique(clusterMemberships1));
if length(B2_rest) == 1
    clusterMemberships2 = 1;
else
    clusterMemberships2 = clusterdata(B2_rest', clustersNo);
end
% clusterMemberships1 = kmeans( A2_rest, clustersNo, 'replicates', 5);
% clusterMemberships2 = kmeans( B2_rest, clustersNo, 'replicates', 5);
clustersNoSet2 = clustersNo;

nextClusterIdx1 = 1;
nextClusterIdx2 = 1;

for j = 1 : clustersNo
    cluster1(j) = clusterMemberships1(nextClusterIdx1);
    cluster2(j) = clusterMemberships2(nextClusterIdx2);
    %clusterMembers1 = find( D_A2 == A2_rest(clusterMemberships1 == cluster1(j)) );
    %clusterMembers2 = find( D_B2 == B2_rest(clusterMemberships2 == cluster2(j)) );
    degreesOfClusterMembers1 = A2_rest(clusterMemberships1 == cluster1(j));
    clusterMembers1 = find( D_A2 <= max(degreesOfClusterMembers1) & ...
        D_A2 >= min(degreesOfClusterMembers1) );
    degreesOfClusterMembers2 = B2_rest(clusterMemberships2 == cluster2(j));
    clusterMembers2 = find( D_B2 <= max(degreesOfClusterMembers2) & ...
        D_B2 >= min(degreesOfClusterMembers2) );
    total_elements = length(clusterMembers1) * length(clusterMembers2);
    Q0(clusterMembers1, clusterMembers2) = 1 ./ ( 1 + relDiffQ(clusterMembers1, clusterMembers2) );
    nextClusterIdx1 = find( clusterMemberships1 == cluster1(j), 1, 'last' )+1;
    nextClusterIdx2 = find( clusterMemberships2 == cluster2(j), 1, 'last' )+1;
    if nextClusterIdx1 > length(clusterMemberships1) || nextClusterIdx2 > length(clusterMemberships2)
        break;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% sort unique degrees of set 1 in descending order

% find the nodes that have not been matched already and operate only on
% them.
total_nodesA1 =  1:1:size(A,2) ;
unmatched_nodesA1 = total_nodesA1(~ismember(total_nodesA1, idx_A1_matched));
total_nodesB1 =  1:1:size(B,2) ;
unmatched_nodesB1 = total_nodesB1(~ismember(total_nodesB1, idx_B1_matched));


A1_deg_unique = sort(unique(D_A1(unmatched_nodesA1)), 'descend');
B1_deg_unique = sort(unique(D_B1(unmatched_nodesB1)), 'descend');
l_A1 = length(A1_deg_unique);
l_B1 = length(B1_deg_unique);

%% match the high degree nodes (up to the knee in the degrees plot)
knee1 = detectKneeAtScreeLikePlot( A1_deg_unique );
knee2 = detectKneeAtScreeLikePlot( B1_deg_unique );
threshold = min( knee1, knee2 );
if threshold == 1
    threshold = 0;
end
for i = 1 : threshold
    idx_A1 = find(D_A1 == A1_deg_unique(i));
    idx_B1 = find(D_B1 == B1_deg_unique(i));
    idx_A1unmatched = ~ismember(idx_A1, idx_A1_matched);
    idx_B1unmatched = ~ismember(idx_B1, idx_B1_matched);
    total_elements = length(idx_A1_unmatched) * length(idx_B1_unmatched);
    
    P0(idx_B1_unmatched, idx_A1_unmatched) = 1 / total_elements;
end

%% cluster the rest of the degrees in clustersNo and initialize uniformly
%% each cluster
A1_rest = A1_deg_unique(threshold+1:l_A1);
B1_rest = B1_deg_unique(threshold+1:l_B1);

%% From the rest nodes do coupled initialization based on the matchings inferred
%% for set 2 (Q).



% lA1_rest = length(A1_rest);
% lB1_rest = length(B1_rest);
% minlength = min(lA1_rest, lB1_rest);
% clustersNo = min( clustersNo, minlength );
% clustersNoSet1 = clustersNo;
% clusterMemberships1 = kmeans( A1_rest, clustersNo, 'replicates', 5);

clustersNoSet1 = 0;
if ( ~isempty(A1_rest) && ~isempty(B1_rest) )
    if length(A1_rest') == 1
        clusterMemberships1 = 1;
    else
        clusterMemberships1 = clusterdata( A1_rest', 1);
    end
    clustersNo = length(unique(clusterMemberships1));
    % clusterMemberships2 = kmeans( B1_rest, clustersNo, 'replicates', 5);
    if length(B1_rest) == 1
        clusterMemberships2 = 1;
    else
        clusterMemberships2 = clusterdata( B1_rest, clustersNo);
    end
    clustersNoSet1 = clustersNo;
    
    nextClusterIdx1 = 1;
    nextClusterIdx2 = 1;
    
    for j = 1 : clustersNo
        cluster1(j) = clusterMemberships1(nextClusterIdx1);
        cluster2(j) = clusterMemberships2(nextClusterIdx2);
        %clusterMembers1 = find( D_A1 == A1_rest(clusterMemberships1 == cluster1(j)) );
        %clusterMembers2 = find( D_B1 == B1_rest(clusterMemberships2 == cluster2(j)) );
        degreesOfClusterMembers1 = A1_rest(clusterMemberships1 == cluster1(j));
        clusterMembers1 = find( D_A1 <= max(degreesOfClusterMembers1) & ...
            D_A1 >= min(degreesOfClusterMembers1) );
        degreesOfClusterMembers2 = B1_rest(clusterMemberships2 == cluster2(j));
        clusterMembers2 = find( D_B1 <= max(degreesOfClusterMembers2) & ...
            D_B1 >= min(degreesOfClusterMembers2) );
        total_elements = length(clusterMembers1) * length(clusterMembers2);
        P0(clusterMembers2, clusterMembers1) = 1 ./ ( 1 + relDiffP(clusterMembers2, clusterMembers1) );
        nextClusterIdx1 = find( clusterMemberships1 == cluster1(j), 1, 'last' )+1;
        nextClusterIdx2 = find( clusterMemberships2 == cluster2(j), 1, 'last' )+1;
        if nextClusterIdx1 > length(clusterMemberships1) || nextClusterIdx2 > length(clusterMemberships2)
            break;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


end