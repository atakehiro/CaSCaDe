function [bkg,rb,rg]=estibkgkmean(I,Nmode)
    % using kmeans classifiction to identify the backgorund
    %
    if nargin==1
        Nmode=2;
    end
    [idx,C]=kmeans(I,Nmode);
    [~,sid]=sort(C,'ascend');
    idx0=zeros(size(idx));
    for k=1:Nmode
        idx0(idx==sid(k))=k;
    end
    Ibg=I(idx0==1);
    rg=[min(I(idx0==1)),max(I(idx0==1))];
    bkg=mean(Ibg);
    rb=std(Ibg);
end