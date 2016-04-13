MAX_ITER=10000;
epsilon=0.000001;
lambda=0.1;
mu=lambda*(size(A,2)^2)/(size(A,1)^2);
cost=zeros(MAX_ITER,1);

%load data
%movie_data=tdfread('ml-100k/u.item','bar');
A=loadMovie('ml-1m/movies.dat');
%try different noise level
noise = [0.01, 0.02, 0.05, 0.01, 0.15, 0.20];
%aggregate structural alike nodes
aggrA=aggr(A);
for level=1:size(noise,2)
    for iter=1:5
        idx=randperm(size(A,1));
        B=A(idx,:);
        B=+xor( B, binornd( ones(size(A,1),size(A,2)), noise(level) ) );
        aggrB=aggr(B);
        tic;
        [permuted_A, final_cost] = GDwINIT(aggrA,aggrB,lambda,mu,MAX_ITER,epsilon);
        time(level,iter)=toc;
        acc(level,iter)=1-sum(sum(abs(aggrB-permuted_A)))/sum(aggrB(:));
    end
end
% or just compute P
P=computeP(A,B,lambda,0.00001);
err=sum(sum(abs(B-P*A)))/(size(A,1)*size(A,2))
% scree plot
rank=zeros(6,1);
for i =1:size(A,1)
   rank(degreeA(i))= rank(degreeA(i))+1;
end

