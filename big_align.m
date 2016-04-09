MAX_ITER=10000;
epsilon=0.000001;
lambda=0.1;
mu=lambda*(size(A,2)^2)/(size(A,1)^2);
cost=zeros(MAX_ITER,1);

%load data
%movie_data=tdfread('ml-100k/u.item','bar');

%try different noise level
noise = [0.01, 0.10, 0.20, 0.30, 0.50, 1.0];
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
for level=1:6
    for iter=1:5
        idx=randperm(size(A,1));
        P_star=eye(size(A,1),size(A,1));
        P_star=P_star(idx,:);
        B=A(idx,:)+noise(level)*randi([0 1], size(A,1),size(A,2));
        tic;
        [permuted_A, final_cost] = GDwINIT(A,B,lambda,mu,MAX_ITER,epsilon);
        time(level,iter)=toc;
        rounded_A=(permuted_A>0.9);
        rounded_B=(B>=1);
        acc(level,iter)=sum(sum(rounded_A==rounded_B))/(size(A,1)*size(A,2));
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

