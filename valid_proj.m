function P=valid_proj(P)
    P(P>1)=1;
    P(P<0)=0;
end