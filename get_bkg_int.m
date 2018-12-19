function [bkg_int]=get_bkg_int(intout0)
    span=100;
    porder=1;
    bkg_int=zeros(size(intout0));
    parfor kbs=1:length(intout0(1,:)) % get mean intensity of each frame,
        bkg_int(:,kbs)=smooth(intout0(:,kbs),span,'sgolay',porder);
    end
end
