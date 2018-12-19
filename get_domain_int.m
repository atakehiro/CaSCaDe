function intout=get_domain_int(im0temp,pattern)
    L=pattern;
    obnum=max(L(:));
    intout=zeros(length(im0temp),obnum); % matrix with size of (frame # * domain #)
    % normalized the intensity shift for correct the bleaching .
    parfor kbs=1:length(im0temp) % get intenisty of each domains at different time frames
        a=(im0temp{kbs});
        for kll=1:obnum
            intout(kbs,kll)=sum(a(L==kll));
        end
    end
    intout=double(intout);
end