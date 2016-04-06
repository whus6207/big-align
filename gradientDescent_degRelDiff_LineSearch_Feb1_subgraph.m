function [ costFrob, cost, accP, accQ, accTotal, accP0, accQ0, accTotal0, zeroElements_P, zeroElements_Q,...
    accPhung, accQhung, accTotalhung, ...
    elapsedTime, mu, lamda, iterations, clustersNoSet1, clustersNoSet2, ...
    distDegA1, distDegA2, distDegB1, distDegB2 ] = ...
                          gradientDescent_degRelDiff_LineSearch_Feb1_subgraph( input1, epsilon, lamda, ...
                                           noiseLevel, input2, clustersNo,...
                                           iter, time, corrSet1, corrSet2, ...
                                           rows, columns, ...
                                           p, outputFile )
%% %%%%%%%%%%%%%%%%%%%%% Gradient Descent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find the matching between two given graphs A, B so that the following  %
%  cost function is minimized:                                            %
%  min_(P,Q) f_{aug} = ||PAQ-B||_F^2 + lamda*sum_{i,j} P_{ij} + ...       %
%                                    + mu*sum_{i,j} Q_{ij}                %
%                                                                         %
%  INPUTS:                                                                %
%      Required: 1) input1: filename of input graph A                     %
%                          (or 'none' if graph will be generated randomly)%
%                          (or 'bundle' if a .mat file with A, B,         %
%                           corrSet1, and corrSet2 will                   %
%                           be loaded. In this case 'input2' is the full  %
%                           path to the .mat file. )                      %
%                2) epsilon: constant that affects convergence rate       %
%                4) lamda: cost for non-zero values in the P matrix       %
%                      mu = (size of P) / (size of Q) * lamda             %
%                5) clustersNo: number of clusters for initialization     %
%                6) iter: max number of iterations                        %
%      Optional: 1) rows: number of rows in the randomly generated        %
%                         bipartite graph A                               %
%                   columns: number of columns in A                       %
%                   p: probability of an edge existing                    %
%                   outputFile: filename to save graph A                  %
%  OUTPUTS: cost of the output matching                                   %
%             +++++                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%etaInit = eta;
changedEntries = 0;
toSkip = 0;
accPhung = 0;
accQhung= 0;
accTotalhung = 0;

%% INPUT READING AND PREPARATION OF THE 2nd GRAPH (if not given)
if (nargin>8)
    [ A, B, corrSet1, corrSet2, toSkip, changedEntries ] = ...
        readAndPrepareInputGraphs( input1, input2, noiseLevel, ...
        corrSet1, corrSet2, ...
        rows, columns, p, outputFile);
else
    [ A, B, corrSet1, corrSet2 ] = readAndPrepareInputGraphs( input1, input2, noiseLevel);
end

[~, name, ext ] = fileparts(input1);
[outputFolder, name2, ext ] = fileparts(input2);


t = cputime;
clustersNoSet1 = 0;
clustersNoSet2 = 0;
distDegA1 = 0;
distDegA2 = 0;
distDegB1 = 0;
distDegB2 = 0;

%% INITIALIZATION OF CORRESPONDENCE MATRICES
[ P0, Q0, clustersNoSet1, clustersNoSet2, distDegA1, distDegA2, distDegB1, distDegB2] = initDegreeClusteredAndPropToRelDeg_coupledSubgraphs(A, B, clustersNo); %initializeBlocksNo(A, B, blocksNo);
Qinit = Q0;
Pinit = P0;
P=P0;
Q=Q0;

% Not both mu and lamda are arbitrary - the user sets lamda, and mu is
% proportional to the ration (size of Q)/(size of P)
ratio = (size(A,2)*size(B,2)) / (size(A,1)*size(B,1));
mu = ratio*lamda;


%% ======================== Sanity Checks =============================

[idxQ, valQ, accP0, accQ0, accTotal0, wrongP, wrongQ ] = ...
                                        verify_subgraphs(P0, Q0, corrSet1, corrSet2);
                                    
accPv(1) = 0;
accQv(1) = 0;
accvTotal(1) = 0;
accPv(2) = accP0;
accQv(2) = accQ0;
accvTotal(2) = accTotal0;

sizeP = size(P,1)*size(P,2);
sizeQ = size(Q,1)*size(Q,2);
zerosP_0iter = sizeP - nnz(P0);
zerosQ_0iter = sizeQ - nnz(Q0);
[costFrob0 cost0] = computeCost(A,B,P0,Q0,lamda,mu);
cost(1) = cost0;
cost(2) = cost0;
costFrob(1) = costFrob0;
costFrob(2) = costFrob0;
k = 2;
relDiff(2) = 1000;

fid = fopen('errors_Jan24.txt','a');

%% Alternating Gradient Descent Algorithm
while ( (k < iter)  &&  ( relDiff(k) > epsilon ) )
    if k > 2
        if cost(k) > 1.1*cost(k-1)
            k
            costInc = 'cost increasing!!!!!!!!!'
            fprintf(fid,'input1: %s \t epsilon: %f\n',input1, epsilon);
            fprintf(fid,'lamda: %f, noiseLevel: %f, input2: %s, clusters: %d \n', lamda, noiseLevel, input2, clustersNo);
            fprintf(fid,'Iter %d: %s\n\n', k, costInc);
        end
    end
    k = k+1;    
    Pold = P;
    derivP = derivativeP(P, Q);
    eta1(k) = lineSearchP(P,Q,derivP);
    P = P - eta1(k) * derivP;
    P = validProjectionP(P);
    derivQ = derivativeQ(P, Q);
    eta2(k) = lineSearchQ(P,Q,derivQ);
    Q = Q - eta2(k) * derivQ;
    Q = validProjectionQ(Q);
    [ costFrob(k) cost(k) ] = computeCost(A,B,P,Q,lamda,mu);
    % VERIFICATION AND ACCURACY COMPUTATION
    [~, ~, accP, accQ, accTotal, wrongP, wrongQ ] = verify_subgraphs(P, Q, corrSet1, corrSet2);
    accPv(k) = accP;
    accQv(k) = accQ;
    accvTotal(k) = accTotal;
    zerosP(k) =  sizeP - nnz(P);
    zerosQ(k) =  sizeQ - nnz(Q);
    relDiff(k) = abs(cost(k-1)-cost(k)); %/cost(k-1);
end

iterations = k;

%% END OF ALGORITHM

elapsedTime = cputime-t

fclose(fid)

%% STATISTICS COMPUTATION
totalElements_P = size(P,1)*size(P,2);
zeroElements_P = totalElements_P - nnz(P);
totalElements_Q = size(Q,1)*size(Q,2);
zeroElements_Q = totalElements_Q - nnz(Q);
percentZeroP = zeroElements_P / totalElements_P;
percentZeroQ = zeroElements_Q / totalElements_Q;
totalEdges = nnz(A);
changedPercent = changedEntries / totalEdges;

output = sprintf('RES_Jan29_LinSearch/RES_LinSearch_cluster0_RelDiff_Feb1_unique_%s_percentNodes_%d_lamda_%d_iter_%d.mat', ...
    name2, noiseLevel, lamda*1000, iterations );

% plot_Cost_AccuracySimple( name2, iterations, cost, costFrob, costToFind, ...
%    costFrobToFind, accPv, accQv, accvTotal, ...
%    [zerosP_0iter,zerosP]/totalElements_P, ...
%    [zerosQ_0iter, zerosQ ]/totalElements_Q, mu, lamda, noiseLevel,...
%    outputFolder, time)%, ...


%% SAVING OUTPUT & PLOTTING COST 'n ACCURACY
save(output , ...
    'totalElements_P', 'zeroElements_P', ...
    'totalElements_Q', 'zeroElements_Q', ...
    'accP', 'accQ', 'eta1', 'eta2',...
    'idxQ', 'valQ', ...
    'wrongP', 'wrongQ', ...
    'P', 'Q', 'iterations', 'elapsedTime', 'totalEdges', ...
    'percentZeroP', 'percentZeroQ', 'Pinit', 'Qinit',...
    'corrSet1', 'corrSet2', 'accPv', 'accQv', 'accvTotal', 'relDiff', ...
    'cost', 'costFrob', 'derivP', 'derivQ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS needed by the gradient descent algorithm

    function [P] = validProjectionP(P)
        P(P>1) = 1;
        P(P<0) = 0;
    end

    function [Q] = validProjectionQ(Q)
        Q(Q>1) = 1;
        Q(Q<0) = 0;
    end

    function [der] = derivativeP( P, Q )
        der = 2*(P*A*Q-B)*Q'*A' + lamda;
    end

    function [der] = derivativeQ( P, Q )
        der = 2*A'*P'*(P*A*Q-B) + mu;
    end

    function [etaNew] = lineSearchP(P,Q,deriv)
        
        daq = deriv*A*Q;
        paq = P*A*Q;
        
        etaNew = ( -2*trace(daq*B') + 2*trace(paq*daq') + lamda*sum(sum(deriv))) / (2*trace(daq*daq'));
    end

    function [etaNew] = lineSearchQ(P,Q,deriv)   
        pad = P*A*deriv;
        paq = P*A*Q;
        
        etaNew = ( 2*trace(paq*pad') - 2*trace(pad*B') + mu*sum(sum(deriv))) / (2*trace(pad*pad'));
    end


end

