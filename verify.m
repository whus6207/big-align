function [ idxQ, valQ, accP, accQ, accTotal, wrongP, wrongQ, stats_P, stats_Q ] = ...
                                           verify(P, Q, corrSet1, corrSet2)
%% %%%%% Evaluating the accuracy of the graph matching algorithm %%%%%%%%%%
%   = checking the found matchings against the ground truth               %
%   *INPUTS*: correspondence matrices P and Q                             %
%          true correspondences between the nodes of the two input graphs %
%   *OUTPUTS*: accP, accQ: accuracy (percentage of correct matchings) for %
%              P and Q respectively                                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[valP, idxP]=sort(P,1,'descend');
% true correspondence AND found correspondence
correctP = 0;
wrongP(1) = 0;
stats_P = 0;
stats_Q = 0;

i=0;
secondMatchP = zeros(size(valP,1),1);
for j = 1 : size(valP,1)
    indicesP = idxP(valP(:,j)==valP(1,j),j);
    trueIdx = corrSet1(j,2);
    indicesSameTopValue = ( indicesP == corrSet1(j,2) );
    numberOfTopProb = sum( indicesSameTopValue );
    if ( valP(:,j) == valP(1,j)*ones(length(valP(:,j)),1) ) % uniform init
        i=i+1;
        wrongP(i) = j;
        %secondMatchP(j) = trueIdx; %idxP(2,j);
        %valP(2,j) = valP( find(idxP(:,j)==trueIdx),j );
    %( sum( idxP(valP(:,j)==valP(1,j),j) == corrSet1(j,2) ) >= 1 ) ) ...
    elseif ( numberOfTopProb >= 1 && valP(1,j) > 0  ) 
       % correctP = correctP + 1*( 1 - max( 0, sum(idxP(valP(:,j)==valP(1,j),j))-1 )/size(P,1) );
        correctP = correctP + 1*( 1 - max(length(indicesP)-1, 0) / size(P,1) );
       % min(1,2*1/sum( idxP(valP(:,j)==valP(1,j),j) ));
        if ( trueIdx ~= idxP(1,j) )
            secondMatchP(j) = trueIdx;
        else
            secondMatchP(j) = idxP(2,j);
        end
    else
        i=i+1;
        wrongP(i) = j;
        secondMatchP(j) = trueIdx; %idxP(2,j);
        %valP(2,j) = valP( find(idxP(:,j)==trueIdx),j );
    end
end
%stats_P = [corrSet1 idxP(1,:)' secondMatchP valP(1,:)' valP(2,:)' ((valP(1,:)-valP(2,:))./valP(1,:))'*100];
accP = correctP / size(valP,1); %size(stats_P,1);

[valQ, idxQ]=sort(Q,2,'descend');
correctQ = 0;
wrongQ(1)=0;

secondMatchQ = zeros(size(valQ,1),1);
i=0;
for j = 1 : size(valQ,1)
    indicesQ = idxQ(j,valQ(j,:)==valQ(j,1));
    trueIdx = corrSet2(j,2);
    indicesSameTopValue = ( indicesQ == corrSet2(j,2) );
    numberOfTopProb = sum( indicesSameTopValue );
    if ( valQ(j,:) == valQ(j,1)*ones(1,length(valQ(j,:))) ) % uniform init
        i=i+1;
        wrongQ(i) = j;
        secondMatchQ(j) = trueIdx;
        %valQ(j,2) = valQ( j, find(idxQ(j,:)==trueIdx) );
   % elseif ( (sum( idxQ(j,valQ(j,:)==valQ(j,1)) == corrSet2(j,2) ) >= 1 ))  ...
   %         && ((valQ(j,1) > 0) )
    elseif ( numberOfTopProb >= 1 && valQ(j,1) > 0 )
        correctQ = correctQ + 1*( 1 - max(length(indicesQ)-1, 0) / size(Q,2) );
        if ( trueIdx ~= idxQ(j,1) )
            secondMatchQ(j) = trueIdx;
        else
            secondMatchQ(j) = idxQ(j,2);
        end
    else
        i=i+1;
        wrongQ(i) = j;
        %secondMatchQ(j) = trueIdx;
        %valQ(j,2) = valQ( j, find(idxQ(j,:)==trueIdx) );
    end
end
%stats_Q = [corrSet2 idxQ(:,1) secondMatchQ valQ(:,1) valQ(:,2) ((valQ(:,1)-valQ(:,2))./valQ(:,1))*100];
%accQ = sum(stats_Q(:,2)-stats_Q(:,3)==0) / size(stats_Q,1);
accQ = correctQ / size(valQ,1); %size(stats_Q,1);

accTotal = ( correctP + correctQ ) / ( size(stats_P,1) + size(stats_Q,1) );

end

