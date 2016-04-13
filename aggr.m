function aggrA=aggr(A)
    i=1;
    aggrA=A;
    while i<=size(aggrA,1)-1
        j=i+1;
        count=1;
        while j<=size(aggrA,1)
           if aggrA(i,:)==aggrA(j,:)
               count=count+1;
               aggrA(j,:)=[];
           else
               j=j+1;
           end
        end
        aggrA(i,:)=aggrA(i,:)*count;
        i=i+1;
    end
end