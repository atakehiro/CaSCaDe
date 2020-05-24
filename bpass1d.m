function res = bpass1d(arr,lb,hb)
    %
    b = double(lb);
    r = round(hb);
    w = 2*r + 1;
    r = ((0:w-1) - r)/(2 * b);
    gx = exp( -r.^2) / (2 * b * sqrt(pi));
    mx = ones(1,w)/w;
    res = arr(:)';
    g = conv(res,gx,'valid');
    tmpres = res;
    res = conv2(tmpres,mx,'valid');
    res0=zeros(size(arr));
    g0 = zeros(size(arr));
    res0((hb+1):end-hb) = res;  % !! wyz
    g0((hb+1):end-hb) = g;
    %res0((lobject+1):end-lobject) = res;
    %g0((lobject+1):end-lobject) = g;
    res=max(g0-res0,0);
end
