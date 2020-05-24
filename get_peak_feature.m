function fs=get_peak_feature(itest)
    
    if length(itest)<5
        itest = [0;0;itest;0;0];
    end
    
    % input fregament I;
    % output features --> fs;
    itest2=itest>mean(itest); % use just sequence
    itest3=(itest(2:end)-itest(1:end-1))>0;
    
    [fr,gof]=fit([1:length(itest)]',itest,'poly1');

    % autocorr
%     try
    [ACF]=autocorr(itest);
    [ACF2]=autocorr(double(itest2));
    [ACF3]=autocorr(double(itest3));
%     catch
%         keyboard
%     end
    
    if isnan(ACF3)
        ACF3=0;
    end
    ACF=[ACF', 0 0 0 0 0 0 ];
    ACF2=[ACF2', 0 0 0 0 0 0 ];
    ACF3=[ACF3', 0 0 0 0 0 0 ];
    
    fs=[];
    k=1;
    fs(k)= length(itest);k=k+1; % this sould be mean. need to correct it.
    fs(k)= max(itest);k=k+1;
    %         fs(k)= find(itest==max(itest))/length(itest);k=k+1;
    fs(k)= min(itest); k=k+1 ;
    fs(k)= fr.p1; k=k+1;
    fs(k)= fr.p2; k=k+1;
    
    fs(k)= gof.rsquare;   k=k+1;
    fs(k)= max(itest)-min(itest) ; k=k+1;
    fs(k)= std(itest) ; k=k+1;
    fs(k)= fs(1)/std(itest) ; k=k+1;
    fs(k)= fs(1)/mean(itest); k=k+1;
    
    fs(k)= skewness(itest); k=k+1;
    fs(k)= kurtosis(itest); k=k+1;
    fs(k)= entropy(itest); k=k+1;
    fs(k)= ACF(2) ; k=k+1;
    fs(k)= ACF(4) ; k=k+1;
    %         fs(k)= ACF0(6) ; k=k+1;
    
    fs(k)= max(0,(ACF(2)-ACF(3))) ; k=k+1;
    fs(k)= max(0,(ACF(3)-ACF(4))) ; k=k+1;
    fs(k)= ACF2(2) ; k=k+1;
    fs(k)= ACF2(4) ; k=k+1;
    fs(k)= max(0,(ACF2(2)-ACF2(3))) ; k=k+1;
    
    fs(k)= max(0,(ACF2(3)-ACF2(4))) ; k=k+1;
    fs(k)= ACF3(2) ; k=k+1;
    fs(k)= ACF3(4) ; k=k+1;
    fs(k)= max(0,(ACF3(2)-ACF3(3))) ; k=k+1;
    fs(k)= max(0,(ACF3(3)-ACF3(4))) ; k=k+1;
    
    
    fs=fs(:);
    
end


