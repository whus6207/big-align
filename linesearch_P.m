function eta_p=linesearch_P(A,B,P,Q,df_dP,lambda)
    dAQ=df_dP*A*Q;
    eta_p=(2*trace((P*A*Q)*dAQ'-dAQ*B')+lambda*sum(df_dP(:)) )/(2*norm(dAQ,'fro')^2);
end