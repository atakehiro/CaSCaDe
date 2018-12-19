function res = bpassW(arr,lnoise,lobject)
    %
    % bandpass filter.
    %
    % Written by Pei-Hsun Wu,
    % Post-doc associate, @ JHU, IMBT.
    %
    b = double(lnoise);
    w = round(lobject);
    N=2*w+1;
    hg=fspecial('gaussian',N, b*sqrt(2));
    ha = fspecial('average',N);
    arra = imfilter(arr,hg-ha,'symmetric','conv');
    rest = max(arra,0);
    res=rest;
end