MAX_ITER=10000;
epsilon=0.000001;
lambda=0.1;
mu=lambda*(size(A,2)^2)/(size(A,1)^2);
cost=zeros(MAX_ITER,1);

%load data
%movie_data=tdfread('ml-100k/u.item','bar');

%aggregate structural alike nodes
i=1;
while i<=size(A,1)-1
    j=i+1;
    while j<=size(A,1)
       if A(i,:)==A(j,:)
           A(j,:)=[];
       else
           j=j+1;
       end
    end
    i=i+1;
end
idx=randperm(size(A,1));
P_star=eye(size(A,1),size(A,1));
P_star=P_star(idx,:);
B=A(idx,:);
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
% or just compute P
P=computeP(A,B,lambda,0.00001);
err=sum(sum(abs(B-P*A)))/(size(A,1)*size(A,2))
% scree plot
rank=zeros(6,1);
for i =1:size(A,1)
   rank(degreeA(i))= rank(degreeA(i))+1;
end