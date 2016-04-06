function [cost,df_dP]=f_aug(A,B,P,Q,lambda,mu)
    PAQ=P*A*Q;
    cost=trace(PAQ*PAQ'-2*PAQ*B'+B*B')+lambda*sum(P(:))+mu*sum(Q(:));
    df_dP=2*(PAQ-B)*Q'*A'+lambda;
end