function [pkout,intoutbw]=peak_detect_v2(intout0,p)
    % pkout : N*4 matrix, (N: number of peaks;
    % 4 columns are 1. domain id, 2.time id, event width, event peak
    % int.
    intoutbw=zeros(size(intout0));
    obnum=size(intout0,2);
    pkout=[];
    for kii=1:obnum
        ytemp=intout0(:,kii);
        bw=ytemp> p.min_int_ed;
        L0=bwlabel(bw);
        bwnew=bw<0; %(zero logical matrix)
        [~,locs]=findpeaks(ytemp,'minpeakheight',p.peak_int_ed,'MINPEAKDISTANCE',p.min_peak_dist_ed);
        for kll=1:length(locs)
            bwnew= bwnew | L0==L0(locs(kll));
        end
        intoutbw(:,kii)=(bwnew(:)); % eliminate signal with less than 5 span away;
        bw=bwnew;
        L0=bwlabel(bw);
        % peak segmentation
        % segment when there are more than one peaks in an detected evnt (bw).
        test=L0(locs);
        kold=1;
        for ktt=1:length(test)-1
            kcurrent=ktt+1;
            if test(kold)==test(kcurrent) % if repeated. then segmented
                id=locs(kold):locs(kcurrent);
                col=find(ytemp(id)==min(ytemp(id)));
                col=col(1);
                h1n=ytemp(id(1));
                h2n=ytemp(id(end));
                h1=ytemp(id(1))-min(ytemp(id));
                h2=ytemp(id(end))-min(ytemp(id));
                hr1=h1/ytemp(id(1));
                hr2=h2/ytemp(id(end));
                d1=col-1;
                d2=length(id)-col;
                thp=3;
                tha=3;
                ccid=col+locs(kold)-1;
                if h1 > thp && h2 > tha && d1>1 && d2>1 && (hr1 >0.5 && hr2 > 0.5) % condiitons for segments.
                    if h1n>p.peak_int_ed*2 && h2n > p.peak_int_ed*2
                        bw(ccid)=0; % apply segment
                        kold=kcurrent; % update the kcurrent value
                    end
                elseif h2>h1
                    kold=kcurrent;
                end
            else
                kold=kcurrent; % update the kcurrent value
            end
        end
        % reanalyze the peak info based on new segmentation result
        L0 = bwlabel(bw); % new label on updated bw;
        % remove the event start from first frame or last til last
        % frame
        temp2=L0(1);
        if temp2>0
            L0(L0==temp2)=0;
        end
        temp2=L0(end);
        if temp2>0
            L0(L0==temp2)=0;
        end
        bw=L0>0; % update binary activation signal
        L0= bwlabel(bw); % update label;
        bwnew = zeros(size(bw));
        % only collect the ones with signal
        for kll=1:max(L0(:))
            % recalculate peak info.
            bwtemp1=L0==kll;
            loctemp=find(bwtemp1==1); % get peaklocation;
            [pkh,ploc]=max(ytemp.*bwtemp1); % peak height and peak locations
            ploc=ploc(1); % if there are more than one peaks with same maximum value, take first one.
            if sum(bwtemp1)> p.min_peak_length
                % only harvest the peak with length more the threshold
                % value.
                pkl=sum(bwtemp1); % peak length
                bwnew= bwnew | bwtemp1; % update binary info
                f_start=loctemp(1);
                f_end=loctemp(end);
                pkres=[kii,ploc,pkl,pkh, f_start f_end];
                pkout=[pkout;pkres];
            end
        end
        intoutbw(:,kii)=(bwnew(:)); % eliminate signal with less than 5 span away;
    end
end