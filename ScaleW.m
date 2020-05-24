function [scaled, Lower, Upper, MaxV, MinV ] = ScaleW(Data, Lower, Upper, MaxV, MinV)
    if (nargin<3)
        Lower = -1;
        Upper = 1;
    elseif (Lower > Upper)
        disp ('Wrong Lower or Upper values!');
    end
    if nargin<=3 % calulate MaxV and MinV
        [MaxV, ~]=max(Data);
        [MinV, ~]=min(Data);
    end
    [R,C]= size(Data);
    scaled=(Data-ones(R,1)*MinV).*(ones(R,1)*((Upper-Lower)*ones(1,C)./(MaxV-MinV)))+Lower;
    scaled=min(Upper, scaled);
    scaled=max(Lower, scaled);
end