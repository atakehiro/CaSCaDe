function res=Cal_anl_main2sa_forreview_x(dat,p)
    % input :
    % im0 : data
    % output :
    % res: resutls output
    %
    % parameter setting;
    %     % basic analysis criteria setting
    %     p.foffset=10; % how many initial frames to exclude in analysis
    %     p.norm_signal='std'; % ('std','bkg','sub') % different way to normalize intenisty
    %     p.spf=1 ; % frame rate at acquisition
    %     % event detection
    %     p.min_int_ed=0.5; % minimum intenisty value for start-end of a event;
    %     p.peak_int_ed=2.0; % minimum peak intesnity value for being considered as signal
    %     p.min_peak_dist_ed=2;
    %     p.min_peak_length=2;
    %     % background trending correction
    %     p.int_correct= 0; % if 1, correct bkg, if 0, no correction.
    %
    % Modified by Yizhi

    T = size(dat,3);    
    im0 = cell(T,1);
    for tt=1:T
        im0{tt} = dat(:,:,tt);
    end
    im0=im0(:);

    % read information regarding images and conditions.
    % system parameter
    knn=1;
    fii=1;
    [hh,ww]=size(im0{1});
    iminfo.Height=hh;
    iminfo.Width=ww;
    
    % acquire image data for this condition
    % set frame of interests for each conditions
    fini=1+p.foffset;
    fend=length(im0);
    frameset0=(fini:fend);
    
    % spatial temporal convolution
    im0temp=im0(frameset0);    
    
    % identify domain candidates
    im3f=bpass3d_v1(im0temp,p);
    [bff]=sum(im3f,3);
    bff=bff/length(frameset0);
    bw=domain_segment(bff,[],p.hb) ;
    L=bwlabel(bw);
    % obtain domains information
    stats=regionprops(bw,'area'); % get area of each mode
    A=[stats.Area];
    obnum=max(L(:)); % all detected node from spatial descrepency

    
    % recording the intenisty of each modes and normalization
    % ---------------------------------------------------------
    
    % initiate output variables
    if obnum>0
        intout=get_domain_int(im0temp,L);
        
        % get normalized intensity by doamin size.(inout0)
        intout0=intout./(ones(length(intout(:,1)),1)*A(:)');
        if p.int_correct % correction for intensity shift, possible photobleaching effect
            [bkg_int]=get_bkg_int(intout0);
            intout0=intout0-bkg_int;
        end
        
        % need to normalize the intout / normalized the intensity.
        bg00=zeros(size(intout0(1,:)))'; rbg00=bg00;
        parfor k99=1:length(intout0(1,:))
            [bg,rbg]=estibkg(intout0(:,k99),7); % estimate the background;
            bg00(k99)=bg;
            rbg00(k99)=rbg;
        end
        
        medmat=ones(size(intout0(:,1)))*bg00';
        stdmat=ones(size(intout0(:,1)))*rbg00';
        % intenity profile with different normalization
        intoutb1=(intout0-medmat)./stdmat; %normalized signal
        intoutb2=(intout0-medmat)./medmat; %normalized signal
        intoutb3=(intout0-medmat); %normalized signal
        switch(p.norm_signal)
            case('std')
                intoutf=intoutb1;
            case('bkg')
                intoutf=intoutb2;
            case('sub')
                intoutf=intoutb3;
        end
        
        % peak detecting at individaul domains and processing.
        [pkout0, intoutbw]=peak_detect_v2(intoutf,p);
        
        % get intenity value of event peak at differnet normalization method
        pk_int1=get_peak_intprofile(pkout0,intoutb1);
        pk_int2=get_peak_intprofile(pkout0,intoutb2);
        pk_int3=get_peak_intprofile(pkout0,intoutb3);
        pkout=[pkout0,[pk_int1 pk_int2 pk_int3]];
        % pkout format: domain id/ peak location/peak length/peak height/
        % peak initial frame/ end frame
    end
    
    
    % get features of peaks
    % ------------------------
    
    %pknum=0;
    if ~isempty(pkout)
        pk_DA=A(pkout(:,1)); % corresponding domain size of individual events
        pk_bg=bg00(pkout(:,1)); % background intensity
        pk_rbg=rbg00(pkout(:,1)); % standard dev.
        pknum=length(pkout(:,1));
        peakfeatures=zeros(pknum,75);
        peak_int_t=cell(pknum,1); % variable to store all detected peak profiles
        
        % go through individual events,
        for kpk=1:pknum
            if mod(kpk,10)==0
                fprintf('Fts %d\n',kpk);
            end
            
            % get peak features for peaks
            intp= intoutf(:,pkout(kpk,1)) ; % intensity profile of this event.
            intp2=intoutbw(:,pkout(kpk,1));
            i1b = bpass1d(intp,1,11);
            i1b2= bpass1d(intp,1,21);
            L1=bwlabel(intp2);
            bwint= L1==L1(pkout(kpk,2));
            itest=intp(bwint);
            itestb=i1b(bwint);
            itestb2=i1b2(bwint);
            fs1=get_peak_feature(itest);
            fs2=get_peak_feature(itestb);
            fs3=get_peak_feature(itestb2);
            peakfeatures(kpk,:)=[fs1',fs2',fs3'];
            peak_int_t{kpk}=itest;
        end
        pkinfo=[knn*ones(pknum,1),fii*ones(pknum,1),1*ones(pknum,1),pkout, pk_DA(:),pk_bg(:), pk_rbg(:)];
    else        
        pkinfo = [];
    end
    
    
    % output
    % ------------
    
    % save organized output resutls
    p.iminfo=iminfo;
    res.param=p; % used parameters setting
    res.frameset0=frameset0; % frames that is used in analysis
    res.bff=bff; % 2D projection image after 3D convolution.
    res.bw=bw; % domain segmetnation image
    res.L=L; % labeled domain segmentation image
    res.A=A; % Area of each domain from domain segmetnation image
    res.obnum=obnum; % number of domain based on domain segmentaion image
    res.medmat=medmat; % median intenisty for each domain
    res.stdmat=stdmat; % standard deviation in intensity profile for each domain
    res.intout=intout; % differnet type of intenisty profiles of domains
    
    res.intoutf=intoutf;
    res.intoutb1=intoutb1;
    res.intoutb2=intoutb2;
    res.intoutb3=intoutb3;
    res.intout0=intout0;
    res.intoutbw=intoutbw; % binarized domain intensity profiles : 1 means detecting as signal. 0 means no signal.
    res.pkinfo=pkinfo; % detected peak informaiton
    res.peakfeatures=peakfeatures; % features of each detected peaks
    res.peak_int_t=peak_int_t;
        
    % gather features
    pkf = res.peakfeatures;
    res.testdata = prepSVMData(pkf);
    
    % get labels
    HW = size(dat,1)*size(dat,2);
    nEvt = size(pkf,1);
    L = res.L;
    roiLst = label2idx(L);
    pks = res.pkinfo;
    evtx = pks(:,4);
    t0 = pks(:,8);
    t1 = pks(:,9);
    evt = cell(0);
    evt2d = cell(0);
    for jj=1:nEvt
        evt2d{jj} = roiLst{evtx(jj)};
        xx = roiLst{evtx(jj)}+HW*((t0(jj):t1(jj))-1);
        evt{jj} = xx(:);
    end
    res.evt = evt;
    res.evt2d = evt2d;
    res.roi = roiLst;
    res.z = pks(:,7);
    res.svm1_pk_class=ones(nEvt,1); % peak goodness based on svm classication
    
end


