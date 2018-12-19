function pk_int=get_peak_intprofile(pkout0,intoutf)
    % obtain peak and sum_intenisty of a given detected signal
    pk_int=[];
    if ~isempty(pkout0)
        pk_int=zeros(size(pkout0(:,1:2)));
        for kk=1:length(pkout0(:,1))
            xid=pkout0(kk,1);
            yid1=pkout0(kk,5):pkout0(kk,6);
            ytemp=intoutf(yid1,xid);
            sii=sum(ytemp); % sum_intenisty
            pkii=max(ytemp); % peak_intensity
            pk_int(kk,:)=[pkii,sii];
        end
    end
end