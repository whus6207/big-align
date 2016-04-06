function [ accPhung, accQhung, accTotalHung, ...
           accP, accQ, accTotal, ...
           matchingPercentQhungInit, matchingPercentQInit, elapsedTime ] = ...
                  twoStepApproach_Feb17( input1, input2, iter, lamda, ...
                          noiseLevel, numberOfEigenValues, hungarianOrNot )
%% %%%%%%%%%%%%%%% Graph matching implementation %%%%%%%%%%%%%%%%%%%%%%%%%%
%   as presented in the appendix                                          %
%     Two-step approach:        %
%     Spectral Clustering, Graph Matching, and Clique Finding''           %
%                                                                         %
% INPUTS: input 1 -> graph A                                              %
%         input 2 -> graph B                                              %
%         iter -> maximum number of iterations                            %
%         epsilon -> small constant for convergence                       %
%         noiseLevel -> level of noise for the random permutation         %
%         numberOfEigenValues -> how many eigenvalues should be computed  %
%                              for the eigendecomposition of the matrices %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wrongInput = 'Input file should be a txt or mat file... Please try again.';

changedEntries = 0;
matchingPercentP = 0;
matchingPercentPhung = 0;
matchingPercentPInit = 0;
matchingPercentPhungInit = 0;
matching = 0;
elapsedTime = 0;
matchingTrans = 0;
correspondenceSet1 = 0;
correspondenceSet2 = 0;
A = 0;
B = 0;
% bipartite graphs
Abip = 0;
Bbip = 0;
k = 0;
almostZero = 10^-8;

%% INPUT READING AND PREPARATION OF THE 2nd GRAPH (if not given)
[~, name, ext ] = fileparts(input1);


if ( strcmp(input1, 'nmf') == 1 || strcmp(ext, '.txt') == 1 || strcmp(ext, '.mat') == 1 ||  strcmp(input1, 'bundle') == 1 || strcmp(input1, 'uni_bundle') == 1 )
    if ( strcmp(input1, 'nmf') == 1)
        load(input2);
        corrSet = correspondenceSet1;
    elseif ( strcmp(input1, 'uni_bundle') == 1)
        load(input2);
    elseif ( strcmp(input1, 'bundle') == 1)
        load(input2);
        Abip = A;
        Bbip = B;
    else
        [ A n1a ] = bipartiteFileToSymmetricAdjacency(input1);
        if ( (nargin>=2) && (strcmp(input2,'none')==1) )
            if (noiseLevel > 0)
                [corrSet1, corrSet2, Bbipartite, changedEntries] = ...
                    randomNoisyBipartiteGraphPermutation(input1, noiseLevel);
            else
                [corrSet1, corrSet2, Bbipartite] = randomBipartiteGraphPermutation(input1, A);
            end
            % convert bipartite graph to unipartite
            [ B n1b ] = bipartiteToSymmetricAdjacency(Bbipartite);
            corrSet2rev = [ corrSet2(:,1)+n1a corrSet2(:,2)+n1b ];
            corrSet = [corrSet1; corrSet2rev];
        else
            [ B n1b ] = bipartiteFileToSymmetricAdjacency(input2);
            corrSet = zeros( 2, max(size(A,1), size(B,1)) );
        end
    end
    
    A = Abip'*Abip;
    B = Bbip'*Bbip;
    t = cputime;
    
    %% INITIALIZATION OF THE CORRESPONDENCE MATRIX
    [ ~, Q0 ] = NMF_init(A,B);
    Qinit = Q0;
    Q=Q0;
    %% ITERATIVE NMF-BASED ALGORITHM: UPDATING THE CORRESPONDENCE MATRIX
    while (k < iter) %(  any(any(abs(Q-Q0) > epsilon)) && (k < iter) )
        k = k + 1;
        a = (Q'*A*Q*B + (Q'*A*Q*B)') / 2;
        denom = Q*a;
        sqrtMat = sqrt( (A*Q*B) ./ (denom + almostZero) );
        Q = Q .* sqrtMat;
        if mod(k,5000) == 0
            k
        end
    end
    
    iterations = k;
    
    
    %% %%%%%% CONVERSION OF CORRESPONDENCE TO PERMUTATION MATRIX %%%%%%%%%%
    % use hungarian algorithm for bipartite graph matching
    % -------------------------------------------------------------------------
    % **Note**: entry (i,j) of Q has the probability of matching node i to node
    % j, while the hungarian algorithm takes as input the cost matrix.
    % Thus, we need to convert the correspondence matrix to cost matrix:
    %        (cost matrix) = max_value_in(corr matrix) - (corr matrix)
    % Problem 1: NaN values in Q ---->
    % my solution: replacing zero denominator with a small value
    % -------------------------------------------------------------------------
    % * 3 inefficient implementations:
    % [ matching, cost ] = Hungarian(max(max(Q))-Q);
    % [ matching, cost ] = munkres(max(max(Q))-Q);
    % [ matching, cost ] = assignmentallpossible(max(max(Q))-Q);
    % I = eye(length(matching), length(matching));
    % matchingMatrix = I(matching,:);
    % -------------------------------------------------------------------------
    
    
    if ( hungarianOrNot == 1)
        [ matchingV, ~ ] =  assignmentsuboptimal2(1-Q);
        I = eye(length(matchingV), length(matchingV));
        matchingQ = I(matchingV,:);
        
        [ P ] = secondStep_appendix( A, B, matchingQ, lamda, numberOfEigenValues );
        P(P<0) = 0;
        
        [ matchingV, ~ ] =  assignmentsuboptimal2(1-P);
        I = speye(length(matchingV), length(matchingV));
        matchingP = I(matchingV,:);
    end
        
    
    if ( hungarianOrNot == 1)
        [ matchingVinit, costHunginit ] =  assignmentsuboptimal2(1-Qinit); 
        I2 = eye(length(matchingVinit), length(matchingVinit));
        matchingInit = I2(matchingVinit,:);
    end
    
        
    %% END OF ALGORITHM
    elapsedTime = cputime-t;
    
    
    %% REFINEMENT STAGE??? - p. 168 in paper
    
    %% VERIFICATION AND ACCURACY COMPUTATION
    if ( hungarianOrNot == 1)
        [~, ~, accPhung, accQhung, accTotalHung, ~, ~, ~, ~] = verify(matchingP, matchingQ, correspondenceSet1, correspondenceSet2);
    end
    
    if ( hungarianOrNot == 1)
        [ matchingPercentQhungInit, ~, ~ ] = NMF_verify(matchingInit', correspondenceSet2);
    end
    
    
    % ------------------------------------------------------------------- %
    % **Note**: In this paper, the authors are computing the transpose of
    % the Q matrix that I have in my gradientDescent program (or the one
    % that is found in Umeyama's approach). To verify the correspondence -
    % using my program's verification module, I need to send the transpose
    % of the Q matrix computed by the current implementation.
    % Meaning of (i,j) entry of matrix Q: the i^th node in the first graph
    % (A) is matched to the j^th node in the second graph.
    % Meaning of (i,j) entry of *MY* or *UMEYAMA's* matrix Q: the j^th node
    % in the first graph (A) is matched to the i^th node in the second graph.
    % ------------------------------------------------------------------- %
    %[ matchingPercentQ, ~, ~ ] = NMF_verify(Q, corrSet2);
    [ matchingPercentQInit, ~, ~ ] = NMF_verify(Qinit', correspondenceSet2);
    %[ matchingPercentP, ~, ~ ] = verify_P(P, corrSet1);
    [~, ~, accP, accQ, accTotal, ~, ~, ~, ~] = verify(P, Q, correspondenceSet1, correspondenceSet2);

    
    %% STATISTICS COMPUTATION
    totalElements_P = size(Q,1)*size(Q,2);
    zeroElements_P = totalElements_P - nnz(Q);
    percentZeroP = zeroElements_P / totalElements_P;
    totalEdges = nnz(A);
    changedPercent = changedEntries / totalEdges;
    
    %% not correct input file format - program ending without executing anything
else
    wrongInput
end

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%initialization of the orthonormal matrix Q
    function [ Q, Q0 ] = NMF_init( A, B )
        %% dimensions of correspondence matrices Q
        p1 = size(B,1);
        p2 = size(A,1);
        
        %% First, initialize the "new" correspondence matrices.
        Q = zeros(p1, p2);
        
        if ( strcmp(numberOfEigenValues, 'all') == 0 )
            [V_A, S_A] = eigs(A, numberOfEigenValues);
            [V_B, S_B] = eigs(B, numberOfEigenValues);
        else
            [V_A, S_A] = eig(full(A));
            [V_B, S_B] = eig(full(B));
        end
        V_Aabs = abs(V_A);
        V_Babs = abs(V_B);
        normalized_A = isNormalized(V_Aabs);
        normalized_B = isNormalized(V_Babs);
        V_Aabs_sorted = sortDescendingOrder(V_Aabs, S_A);
        V_Babs_sorted = sortDescendingOrder(V_Babs, S_B);
        % The rows of V_Aabs and V_Babs should be non-negative vectors with length 1.
        Q0 = V_Aabs_sorted * V_Babs_sorted';
    end

    function [ sortedEigV, Sdiag] = sortDescendingOrder(V, S)
        [ Sdiag idx ]  = sort(diag(S),'descend');
        sortedEigV = V(:,idx);
    end

    function [ normalizedMatrix ] = normalizeMatrixPerRow( M )
        % square root of sum of squares of elements per row
        sqrtOfSumOfSquares = diag(sqrt(M*M'));
        % number of rows
        r = size(M,1);
        % number of columns
        c = size(M,2);
        normalizedMatrix = M ./ repmat(sqrtOfSumOfSquares, 1, c);
        correct = 0;
        for i = 1 : r
            if ( norm(normalizedMatrix(i,:),2) > 0.999 && norm(normalizedMatrix(i,:),2) < 1.0001 )
                correct = correct + 1;
            else
                norm(normalizedMatrix(i,:),2)
            end
        end
        correctPercent = correct / r;
    end


    function [ orthonormal ] = isOrthonormal( M )
        unitMat = eye(size(M,1));
        if ( M*M' == unitMat && M'*M == unitMat )
            orthonormal = 1;
        else
            orthonormal = 0;
        end
    end

    function [ normalized ] = isNormalized( M )
        max_dimension = max(size(M,1), size(M,2));
        normalized = ones(2, max_dimension);
        % checking the rows
        for i = 1 : size(M,1)
            if ( norm(M(i,:),2) < 0.9999 || norm(M(i,:),2) > 1.0001 )
                normalized(1,i) = norm(M(i,:),2);
            end
        end
        % checking the columns
        for j = 1 : size(M,2)
            if ( norm(M(:,j),2) < 0.9999 || norm(M(:,j),2) > 1.0001 )
                normalized(2,j) = norm(M(:,j),2);
            end
        end
    end

end