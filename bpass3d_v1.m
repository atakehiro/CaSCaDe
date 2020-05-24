function im3f=bpass3d_v1(im,param)
    % 3D band pass convolution
    % im : stacked images
    % format1 : N*1 cells and each cell corrposdence to a image with size
    % of h*w
    % format2 : h*w*N matrices
    % % processing is done using single precision- for reducing memory loading
    % bim : processed image
    % developed by : Pei-Hsun Wu, Ph.D
    % 10/28/2014 @ Johns Hopkins University
    % set bandpass process parameteres
    if nargin==1
        lb=1; % low bound size for in-x,y dim
        hb=11; % high bound size for in-x,y dim
        zlb=1;% low bound size for in-z(t) dim
        zhb=21; % high bound size for in-z(t) dim
    else
        lb=param.lb;
        hb=param.hb;
        zlb=param.zlb;
        zhb=param.zhb;
    end
    % obtain the 3D image size
    if iscell(im)
        znum=length(im);
        [h,w]=size(im{1});
    else
        [h,w,znum]=size(im);
    end
    im3=single(zeros(h,w,znum));
    
    % convolution in 2D (x,y direction first through whole stack
    if iscell(im)
        parfor kbs=1:znum % time period without drug
            a=single(im{kbs});
            b=bpassW(a,lb,hb);
            im3(:,:,kbs)=single(b);
        end
    else
        parfor kbs=1:znum % time period without drug
            a=single(im(:,:,kbs));
            b=bpassW(a,lb,hb);
            im3(:,:,kbs)=single(b);
        end
    end
    im3f=im3;
    % bandpass filtering over z (t) dim
    parfor kx=1:w
    %for kx=1:w
        for ky=1:h
            temp=squeeze(im3(ky,kx,:));
            im3f(ky,kx,:)=single(bpass1d(temp,zlb,zhb));
        end
    end
end




