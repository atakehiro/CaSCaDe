function [bkg,rb]=estibkg(I,iter,gmode)
    % estimate the image background using a iterative process
    %
    %*******************************************************
    % written by :
    % Pei-Hsun Wu, PhD
    % Institue for Nano-bio technology
    % Johns Hopkins University
    %
    % Last update: 08/19/2013
    %********************************************************%
    %
   
    if nargin==2
        gmode='iterNspace';
    end
    drr=2;
    binnum=40;
    switch(gmode)
        case('iter')
            I=double(I(:));
            itop= mean(I)+ drr* std(I) ; % 16 bit
            for mm=1:iter
                [c1,c2]=hist(I(I < itop),40);
                mimg= c2(c1 == max(c1)) ; mimg=mimg(1);
                ira=std(I(I < itop));
                itop= mimg + drr*ira ;
            end
            bkg=mimg ;
            rb= ira ;
        case('gaufit')
            I=double(I(:));
            mimj = mode(I); % initial guess
            rb=std(I);
            [c1,c2]=hist((I((I < mimj +3*rb))),30);
            [fr]=fit(c2(:),c1(:),'gauss1');
            bkg=fr.b1 ;
            rb=fr.c1/sqrt(2) ;
        case('dir')
            I=double(I(:));
            bkg= mean(I);
            rb=std(I);
        case('iterNspace')
            I1=I;
            I=double(I(:));
            
            itop= mean(I)+ drr* std(I) ; % 16 bit
            for mm=1:iter
                [c1,c2]=hist(I(I < itop),binnum);
                mimg= c2(c1 == max(c1)) ; mimg=mimg(1);
                ira=std(I(I < itop));
                itop= mimg + drr*ira;
            end
            bkg=mimg;
            cc=I1 < bkg;
            cc=imclose(cc,strel('disk',3));
            bkg= mean(I1(cc));
            rb=std(I1(cc));
    end
    
    % !! add by Yizhi. The above noise estimation may sometimes fail
    dx = (sqrt((I(2:end) - I(1:end-1)).^2/2))/0.6743;
    rb = median(dx(:));
    
end



