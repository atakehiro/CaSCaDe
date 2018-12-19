function bw=domain_segment(bff,Areacut,hb)
    % segment the calcium domains area
    if nargin==1 || isempty(Areacut)
        Areacut=25;
    end
    [hh,ww]=size(bff);
    if hh==ww
        mskc=mskcircle(length(bff));
    else
        rr=5;
        mskc=zeros(size(bff));
        mskc(rr+1:end-rr,rr+1:end-rr)=1;
    end
    test=bff.*mskc;
    [bge,rbge]=estibkgkmean(test(test>0));
    bw=bff>bge+2*rbge;
    % bw=bff>80; %(100) set threshold to binarize the locations demonstrate spatial temporal distintive
    % xyc is the just to highlighted the center region of images. morelikely where the cell is located
    %b000=bpassW(bff,3,61); % default 3,21;
    b000=bpassW(bff,3,hb); % default 3,21;
    
    nbw=imregionalmax(b000,8);
    nbw=imdilate(nbw,strel('disk',2));
    nbw=nbw & bw ;
    % watershed segmentation
    cbw=bw | nbw;
    D=bwdist(nbw);
    DL=watershed(D);
    cbw(DL==0)=0;
    cbw00=cbw<0; % preallocateing cbw00
    ln=bwlabel(nbw);
    lc=bwlabel(cbw);
    DLrev=DL < 0 ;
    for k=1:max(lc(:))
        idx=ln(lc==k);
        [uidx]=unique(idx);
        if length(uidx)==1
            bwtemp=lc==k;
            cbw00=cbw00 | bwtemp; % collect these non-nucleated reguion
            bwtemp=imdilate(bwtemp,strel('disk',2));
            DL00= DL==0 & imdilate(bwtemp,strel('disk',2));
            DLrev=DL00 | DLrev;
        end
    end
    bw=cbw & ~cbw00; % does not take the one without nuclesate
    bw=bwareaopen(bw,Areacut); % get rid of nodes with small area size
    bw=imfill(bw,'hole');
end