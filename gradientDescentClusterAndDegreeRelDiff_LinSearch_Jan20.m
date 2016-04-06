function [ costFrob, cost, accP, accQ, accTotal, accP0, accQ0, accTotal0, zeroElements_P, zeroElements_Q,...
    accPhung, accQhung, accTotalhung, ...
    elapsedTime, mu, lamda, iterations, clustersNoSet1, clustersNoSet2, ...
    distDegA1, distDegA2, distDegB1, distDegB2 ] = ...
                          gradientDescentClusterAndDegreeRelDiff_LinSearch_Jan20( input1, epsilon, lamda, ...
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
%                3) eta: step for gradient descent                        %
%                4) lamda: cost for non-zero values in the P matrix       %
%                      mu = (size of P) / (size of Q) * lamda             %
%                5) blocksNo: number of blocks for initialization         %
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
%[P, Q, P0, Q0, clustersNoSet1, clustersNoSet2, distDegA1, distDegA2, distDegB1, distDegB2] = ...
% [P, Q, P0, Q0] = initializeBlocksNo(A, B, clustersNo); % clustersNo = blocksNo
%[P, Q, P0uncoupled, Q0uncoupled, clustersNoSet1, clustersNoSet2, distDegA1, distDegA2, distDegB1, distDegB2] = initDegreeClustered(A, B, clustersNo); %initializeBlocksNo(A, B, blocksNo);
%[~, ~, P0diff, Q0diff, clustersNoSet1, clustersNoSet2, distDegA1, distDegA2, distDegB1, distDegB2] = initDegreeClustered_coupled(A, B, clustersNo);
[~, ~, P0, Q0, clustersNoSet1, clustersNoSet2, distDegA1, distDegA2, distDegB1, distDegB2] = initDegreeClusteredAndPropToRelDeg_coupled(A, B, clustersNo); %initializeBlocksNo(A, B, blocksNo);
Qinit = Q0;
Pinit = P0;
P=P0;
Q=Q0;
%[P, Q, P0, Q0] = initialize(A,B);
%[P, Q, P0, Q0] = initializeMergeDegrees(A, B, blocksNo);

% Not both mu and lamda are arbitrary - the user sets lamda, and mu is
% proportional to the ration (size of Q)/(size of P)
ratio = (size(A,2)*size(B,2)) / (size(A,1)*size(B,1));
%if (ratio*lamda) > 1
%    mu = min(ratio*lamda, lamda*100);
%else
%    mu = ratio*lamda;
%end
mu = ratio*lamda;


%% ======================== Sanity Checks =============================
% Computing the real cost - the one that we ideally want to find
pdim = length(corrSet1);
Preal = zeros(pdim, pdim);
for i = 1 : pdim
    Preal(corrSet1(i,2), corrSet1(i,1)) = 1;
end
qdim = length(corrSet2);
Qreal = zeros(qdim, qdim);
for i = 1 : qdim
    Qreal(corrSet2(i,1), corrSet2(i,2)) = 1;
end
[costFrobToFind costToFind] = computeCost(A,B,Preal,Qreal,lamda,mu)
% Checking if the derivative is zero
derivP = derivativeP( Preal, Qreal );
derivQ = derivativeQ( Preal, Qreal );

[idxQ, valQ, accP0, accQ0, accTotal0, wrongP, wrongQ, stats_P, stats_Q] = ...
                                        verify(P0, Q0, corrSet1, corrSet2);
                                    
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

fid = fopen('errors_Jan20.txt','a');

%% Alternating Gradient Descent Algorithm
while ( (k < iter)  &&  ( relDiff(k) > epsilon ) )                   %( ( any(any(abs(P-P0) > epsilon)) || any(any(abs(Q-Q0) > epsilon))) && (k < iter) )
    if k > 2
        if cost(k) > cost(k-1)
            costInc = 'cost increasing!!!!!!!!!'
            fprintf(fid,'input1: %s \t epsilon: %f\n',input1, epsilon)
            fprintf(fid,'lamda: %f, noiseLevel: %f, input2: %s, clusters: %d \n', lamda, noiseLevel, input2, clustersNo)
            fprintf(fid,'%s \n\n', costInc);
        end
    end
    k = k+1;    
%     Q = Q - eta * derivativeQ(P, Q);
    Pold = P;
    derivP = derivativeP(P, Q);
    eta1(k) = linSearchP(P,Q,derivP);
    P = P - eta1(k) * derivP;
    P = validProjectionP(P);
    derivQ = derivativeQ(P, Q);
    eta2(k) = linSearchQ(P,Q,derivQ);
    Q = Q - eta2(k) * derivQ;
    Q = validProjectionQ(Q);
    [ costFrob(k) cost(k) ] = computeCost(A,B,P,Q,lamda,mu);
    [~, ~, accP, accQ, accTotal, wrongP, wrongQ, stats_P, stats_Q] = verify(P, Q, corrSet1, corrSet2);
    accPv(k) = accP;
    accQv(k) = accQ;
    accvTotal(k) = accTotal;
    zerosP(k) =  sizeP - nnz(P);
    zerosQ(k) =  sizeQ - nnz(Q);
    %eta = updateEta(k+1);
    relDiff(k) = abs(cost(k-1)-cost(k))/cost(k-1);
end

iterations = k;

derivP = derivativeP( P, Q );
derivQ = derivativeQ( P, Q );

% Cost Computation
%[ costFrob(k) cost(k) ] = computeCost(A,B,P,Q,lamda,mu);
%costThreshold = comparePermutedGraph1AndGraph2(A,B,P,Q,lamda,mu);

%% END OF ALGORITHM
elapsedTime = cputime-t

fclose(fid)

%% VERIFICATION AND ACCURACY COMPUTATION
% Translation of correspondence matrices to node/community matching

% if (toSkip == 0)
%     matching based on maximum probability (our approach)
%     [idxQ, valQ, accP, accQ, accTotal, wrongP, wrongQ, stats_P, stats_Q] = ...
%     verify(P, Q, corrSet1, corrSet2);
% else
%     [valP, idxP] = sort(P,1,'descend');
%     [valQ, idxQ] = sort(Q,2,'descend');
%     stats_P = [corrSet1' idxP(1:2,:)' (valP(1,:)-valP(2,:))'];
%     stats_Q = [corrSet2' idxQ(:,1:2) (valQ(:,1)-valQ(:,2))];
%     accP = 0;
%     accQ = 0;
%     wrongP = 0;
%     wrongQ = 0;
% end

%% STATISTICS COMPUTATION
totalElements_P = size(P,1)*size(P,2);
zeroElements_P = totalElements_P - nnz(P);
totalElements_Q = size(Q,1)*size(Q,2);
zeroElements_Q = totalElements_Q - nnz(Q);
percentZeroP = zeroElements_P / totalElements_P;
percentZeroQ = zeroElements_Q / totalElements_Q;
totalEdges = nnz(A);
changedPercent = changedEntries / totalEdges;

output = sprintf('RES_Jan19_LinSearch/RES_LinSearch_cluster0_RelDiff_Jan20_unique_%s_noise_%d_lamda_%d_iter_%d.mat', ...
    name2, noiseLevel, lamda*1000, iterations );

% plot_Cost_AccuracySimple( name2, iterations, cost, costFrob, costToFind, ...
%    costFrobToFind, accPv, accQv, accvTotal, ...
%    [zerosP_0iter,zerosP]/totalElements_P, ...
%    [zerosQ_0iter, zerosQ ]/totalElements_Q, mu, lamda, eta, noiseLevel,...
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
    'costFrobToFind', 'costToFind', 'cost', 'costFrob', 'derivP', 'derivQ')


%%%accPhung, accQhung, accTotalhung )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTIONS needed by the gradient descent algorithm
%     function [P, Q, P0, Q0] = initializeUniform(A, B)
%         pDim1 = size(B,1);
%         qDim1 = size(A,1);
%         pDim2 = qDim1;
%         qDim2 = pDim1;
%         P0 = ones(pDim1, pDim2) ./ pDim2;
%         Q0 = ones(qDim1, qDim2) ./ qDim2;
%         P = zeros(pDim1, pDim2);
%         Q = zeros(qDim1, qDim2);
%     end

    function [P] = validProjectionP(P)
        P(P>1) = 1;
        P(P<0) = 0;
    end

    function [Q] = validProjectionQ(Q)
        Q(Q>1) = 1;
        Q(Q<0) = 0;
    end

%     function [etaNew] = updateEta( iter )
%         beta = 0.9;
%         etaNew = etaInit*beta^(iter-1);
%         etaNew = eta;
%         if iter > 1000
%            etaNew = eta/sqrt(iter);
%        end
%     end

    function [der] = derivativeP( P, Q )
        der = 2*(P*A*Q-B)*Q'*A' + lamda;
    end

    function [der] = derivativeQ( P, Q )
        der = 2*A'*P'*(P*A*Q-B) + mu;
    end

%     function [val] =  objFunP( P1, aqb, aqqa, sumP, sumQ )
%         val = trace((P1'*P1)*aqqa)-2*trace(P1*aqb)+lamda*sumP+mu*sumQ;
%     end

    function [etaNew] = linSearchP(P,Q,deriv)
        
        daq = deriv*A*Q;
        paq = P*A*Q;
        
        etaNew = ( -2*trace(daq*B') + 2*trace(paq*daq') + lamda*sum(sum(deriv))) / (2*trace(daq*daq'));
%         
%         eta = [0:0.00001:0.1];
%        
%         for j = 1 : length(eta)
%               P1 = P-eta(j)*deriv;
%               objFun(j) = trace((P1)*A*Q*(P1*A*Q)') - 2*trace((P1*A*Q)*B') + lamda*sum(sum(P1)) + mu*sum(sum(Q));
%               %             objFunVal(j) = objFunP(P1, aqb, aqqa, sumP, sumQ);
%         end
%         
%         plot(eta, objFun);
        %etaNew = ( trace(daq*B') + trace(paq*daq') ) / trace(daq*daq');
%         number = 500;
%         etaVal = linspace(0.000001, 0.1, number);
%         aq = A*Q;
%         aqqa = aq*aq';
%         aqb = aq*B';
%         sumQ = sum(sum(Q));
%         for j = 1 : number
%             P1 = P - etaVal(j)*deriv;
%             sumP = sum(sum(P1));
%             objFunVal(j) = objFunP(P1, aqb, aqqa, sumP, sumQ);
%         end
%         [~,idx] = min(objFunVal);
%         etaNew = etaVal(idx);
    end


%     function [val] =  objFunQ( Q1, bpa, appa, sumP, sumQ )
%         val = trace(P*aqqa)-2*trace(P1*aqb)+lamda*sumP+mu*sumQ;
%     end

    function [etaNew] = linSearchQ(P,Q,deriv)
        
        pad = P*A*deriv;
        paq = P*A*Q;
        
        etaNew = ( 2*trace(paq*pad') - 2*trace(pad*B') + mu*sum(sum(deriv))) / (2*trace(pad*pad'));
        
%         eta = [0:0.00001:0.1];
%        
%         for j = 1 : length(eta)
%               Q1 = Q-eta(j)*deriv;
%               objFun(j) = trace(P*A*Q1*(P*A*Q1)') - 2*trace((P*A*Q1)*B') + lamda*sum(sum(P)) + mu*sum(sum(Q1));
%               %             objFunVal(j) = objFunP(P1, aqb, aqqa, sumP, sumQ);
%         end
%         
%         plot(eta, objFun);
    end

% function [val] =  objFunQ(Q1, bpa, appa, sumP, sumQ)
%         val = trace(appa*(Q1*Q1'))-2*trace(bpa*Q1)+lamda*sumP+mu*sumQ;
%         %trace(PAQ(PAQ)'-2PAQB')+lamda*sumP+mu*sumQ;
%     end
% 
%     function [etaNew, idx] = linSearchQ(P,Q,deriv)
%         number = 500;
%         etaVal = linspace(0.000001, 0.1, number);
%         pa = P*A;
%         appa = pa'*pa;
%         bpa = B'*pa;
%         sumP = sum(sum(P));
%         for j = 1 : number
%             Q1 = Q - etaVal(j)*deriv;
%             sumQ = sum(sum(Q1));
%             objFunVal(j) = objFunQ(Q1, bpa, appa, sumP, sumQ);
%         end
%         [~,idx] = min(objFunVal);
%         etaNew = etaVal(idx);
%     end


end

