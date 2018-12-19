function msk=mskcircle(sz)
    m = sz;
    x = 0:(m-1) ;
    cent = (m-1)/2;
    x2 = (x-cent).^2;
    dst=zeros(m,m);
    for i=1:m
        dst(i,:)=sqrt((i-1-cent)^2+x2);
    end
    
    ind=find(dst <= cent);
    
    msk=zeros([m,m]);
    msk(ind)=1.0;
end