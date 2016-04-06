function P=computeP(A,B,lambda,err)
    USU=pinv(A*A',err);
    P=B*A'*USU-lambda/2*ones(size(A,1))*USU;
end