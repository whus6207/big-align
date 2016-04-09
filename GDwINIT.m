function [permuted_A, final_cost] = GDwINIT(A,B,lambda,mu,MAX_ITER,epsilon)
    cost=zeros(MAX_ITER,1);
    %NET INITialize
    degreeA=sum(A,2);
    [sortedDegA, sortedIdA]=sort(degreeA,'descend');
    degreeB=sum(B,2);
    [sortedDegB, sortedIdB]=sort(degreeB,'descend');
    P=zeros(size(A,1),size(A,1));
    for i = 1:size(A,1)
        P(sortedIdB(i),sortedIdA(i))=1;
    end
    Q=eye(size(A,2),size(A,2));
    %alternating projected gradient descent
    k=2;
    [cost(2), df_dP]=f_aug(A,B,P,Q,lambda,mu);
    %cost(2)=computeCost(A,B,P,Q,lambda,mu);
    while abs(cost(k)-cost(k-1))/cost(k-1)>epsilon && k<MAX_ITER
        if k>2 && cost(k)>cost(k-1)
           costInc='cost increasing!!'
        end
        k=k+1;
        eta_p=linesearch_P(A,B,P,Q,df_dP,lambda);
        P=P-eta_p*df_dP;
        P=valid_proj(P);
        [cost(k), df_dP]=f_aug(A,B,P,Q,lambda,mu);
        %[ costFrob(k) cost(k) ]=computeCost(A,B,P,Q,lambda,mu);
    end
    final_cost=cost(k);
    permuted_A=P*A;
end