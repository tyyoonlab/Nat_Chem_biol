disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%  Force-ramp Measurement  %%%%%%%%%%%%%%%%%%%%%%%%%%%']);
a1 = 1;
c1 = 1;
for i =1:100
    pause;
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    disp(['[Click] <-- : Analyzing again previous cycle'])
    disp(['[Click] --> : Analyzing FEC!'])
    disp(['[Click] esc : Stop analyzing'])
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    stateanalysis_seq0 = getkey;
    
    if stateanalysis_seq0 == 28
        if i == 1
            i = 1;
        else
            i = i - 1;
        end
    elseif stateanalysis_seq0 == 29
        i = i;
    else
        return
    end
    disp(['%%%%%%%%%%%%%',num2str(i),'# cycle analysis start!%%%%%%%%%%%%%%%']);
    %% Sectioning each FEC cycle
    dcheck=mean(Force{1}(round(a1):round(c1)));
    figure(99);
    axeshandles = gca;
    hold on;
    plot(a1./mean(fps{1}),dcheck, 'rx', 'linewidth', 3, 'markersize', 20);
    plot(c1./mean(fps{1}),dcheck, 'rx', 'linewidth', 3, 'markersize', 20);

    pause;
    getpts
    tisf{i} = ans;
    tisfconvert{i} = tisf{i}.*mean(fps{1});
    a1 = tisfconvert{i}(1);
    c1 = tisfconvert{i}(end);
    
    FEforce{i} = Force{1}(tisfconvert{i}(1):tisfconvert{i}(end));
    
    if 1
        pause;
        getpts
        tisz{i} = ans;
        tiszconvert{i} = tisz{i}.*mean(fps{1});
        temdfz{i} = find(FEforce{i} == max(FEforce{i}));
        dzinfe(i) = tisfconvert{i}(1)+temdfz{i}(round(length(temdfz{i})/2)) - tiszconvert{i}(end);
    end
    
    FEtracez{i} =Trace{1}(tisfconvert{i}(1)-dzinfe(i):tisfconvert{i}(end)-dzinfe(i),3);
    FEtracez_cor{i} = Tracez_cor{1}(tisfconvert{i}(1)-dzinfe(i):tisfconvert{i}(end)-dzinfe(i));
    FETracex{i} =Trace{1}(tisfconvert{i}(1)-dzinfe(i):tisfconvert{i}(end)-dzinfe(i),1);
    FETracey{i} =Trace{1}(tisfconvert{i}(1)-dzinfe(i):tisfconvert{i}(end)-dzinfe(i),2);
    
   %% Force - linear interpolation
    figure(100+i);hold on;
    subplot(2,1,1)
    plot(1:length(FEforce{i}),FEforce{i},'r','LineWidth',1.5);
    
    subplot(2,1,2); hold on;
    plot(1:length(FEtracez{i}),FEtracez{i},'k','LineWidth',0.5);
    plot(1:length(FEtracez_cor{i}),FEtracez_cor{i},'r','LineWidth',0.5);
    
    disp(['%%%%%%%%%%%%%%% Linear interpolation of force and aligning FECs %%%%%%%%%%%%%%']);
    
    zmins = find(round(FEforce{i}*10)/10 == 1);
    if length(zmins) == 1
        findzmin{i}(1) = zmins;
    else
        findzmin{i} = zmins;
    end
    
    figure(100+i);hold on;
    onef = linspace(1,1,length(findzmin{i}));
    sc=10;
    subplot(2,1,1);hold on;
    scatter(findzmin{i},onef,sc,'bo','filled');
    subplot(2,1,2); hold on;
    plot(findzmin{i},FEtracez{i}(findzmin{i}),'b','LineWidth',1);
    
    tracesync(i) = mean(FEtracez{i}(findzmin{i})) + zdat_model_Nstate(find(Fdat_model == 1));
    
    FEtracez{i} = FEtracez{i} - tracesync(i);
    FEtracez_cor{i} = FEtracez_cor{i} - tracesync(i);
     
    display(['Check force-ramp regime'])
    
    pause;
    if 0
        getpts
        tiszf{i} = ans;
        tisctoh{i} = round(tiszf{i})
        FEforcecorr{i} = FEforce{i};
        
        if FEforce{i}(tisctoh{i}(end)) > FEforce{i}(tisctoh{i}(2))
            stretchingf{i} = linspace(1,round(FEforce{i}(tisctoh{i}(end))*10)/10,length(tisctoh{i}(1):tisctoh{i}(end)));
            FEforcecorr{i}(1:tisctoh{i}(1)-1) = FEforce{i}(1:tisctoh{i}(1)-1);
            FEforcecorr{i}(tisctoh{i}(1):length(tisctoh{i}(1):tisctoh{i}(end))+tisctoh{i}(1)-1) = stretchingf{i};
        else
            stretchingf{i} = linspace(1,round(FEforce{i}(tisctoh{i}(2))*10)/10,length(tisctoh{i}(1):tisctoh{i}(2)));
            % If there is kink force data;
            %relaxingf{i} = FEforce{i}(tisctoh{i}(2):tisctoh{i}(end));
            relaxingf{i} = linspace(round(FEforce{i}(tisctoh{i}(2))*10)/10,FEforce{i}(tisctoh{i}(end)),length(tisctoh{i}(2):tisctoh{i}(end)));
            FEforcecorr{i}(1:tisctoh{i}(1)-1) = FEforce{i}(1:tisctoh{i}(1)-1);
            FEforcecorr{i}(tisctoh{i}(1):length(tisctoh{i}(1):tisctoh{i}(2))+tisctoh{i}(1)-1) = stretchingf{i};
            FEforcecorr{i}(tisctoh{i}(2)+1:length(tisctoh{i}(2):tisctoh{i}(end))+tisctoh{i}(2)) = relaxingf{i};
        end
    end
    FEforcecorr{i} =  smooth(FEforce{i},1200,'moving');
    
    close (100+i);
    figure(100+i);hold on;
    plot(1:length(FEforce{i}),FEforce{i},'k','LineWidth',1);
    plot(1:length(FEforcecorr{i}),FEforcecorr{i},'r','LineWidth',1);
    
    %% Align FEC curves with eWLC model
    
    display(['%%%%%%%%%%%%% Translation of extension for alignment %%%%%%%%%%%%']) 
        
    smoothFEtracez{i} = medfilt1(FEtracez{i},mean(fps{1})/5);
    smoothFEtracez_cor{i} = medfilt1(FEtracez_cor{i},mean(fps{1})/5);
    
    m = 2;
    x = 20;
    f=figure(200+i);
    Screensize = get( groot, 'Screensize' );
    set(f, 'Position', [300, 125, Screensize(3)*0.8, Screensize(4)*0.75 ]);
    figp = get(gcf,'Position');
    set(0,'DefaultFigurePosition', figp);
    axis([-200 250 0 50])
    grid on;
    xlabel('z (nm)'); ylabel('Force(pN)');
    
    while 1
        f; hold on;
        figure(200+i); 
        plot(FEtracez{i},FEforcecorr{i},'color',[0.5 0.5 0.5],'LineWidth',0.1);
        plot(FEtracez_cor{i},FEforcecorr{i},'color',[1 0.5 0.5],'LineWidth',0.1);
        plot(zdat_model_plusMP_pp,Fdat_model,'c','LineWidth',2,'LineStyle','--');
        hold on;
        plot(zdat_model_Nstate,Fdat_model,'r','LineWidth',2,'LineStyle','--');
        hold on;
        plot(zdat_model_plusMP_hs,Fdat_model,'m','LineWidth',2,'LineStyle','--')
        for ij = 1:numel(Iu_ratio)
             plot(zdat_model_plusMP_interU{ij},Fdat_model,'g','LineWidth',1,'LineStyle','--')
        end
        plot(smoothFEtracez{i},FEforcecorr{i},'k','LineWidth',0.5);
        plot(smoothFEtracez_cor{i},FEforcecorr{i},'r','LineWidth',0.5);
        ans = getkey
        if ans == 28
            FEtracez{i} = FEtracez{i} - m;
            smoothFEtracez{i} = smoothFEtracez{i} - m;
            FEtracez_cor{i} = FEtracez_cor{i} - m;
            smoothFEtracez_cor{i} = smoothFEtracez_cor{i} - m;
        elseif ans == 29
            FEtracez{i} = FEtracez{i} + m;
            smoothFEtracez{i} = smoothFEtracez{i} + m;
            FEtracez_cor{i} = FEtracez_cor{i} + m;
            smoothFEtracez_cor{i} = smoothFEtracez_cor{i} + m;
        elseif ans == 30
            FEtracez{i} = FEtracez{i} + x;
            smoothFEtracez{i} = smoothFEtracez{i} + x;
            FEtracez_cor{i} = FEtracez_cor{i} + x;
            smoothFEtracez_cor{i} = smoothFEtracez_cor{i} + x;
        elseif ans == 31
            FEtracez{i} = FEtracez{i} - x;
            smoothFEtracez{i} = smoothFEtracez{i} - x;
            FEtracez_cor{i} = FEtracez_cor{i} - x;
            smoothFEtracez_cor{i} = smoothFEtracez_cor{i} - x;
        else
            break
        end
        clf(f);
    end
    htocc = 0;
    if htocc == 1
        
        [zu,fu] = getline;
        
        uz{i} = zu;
        uf{i} = fu;
        
        A0 = 67.00787;
        lambda0 = 2.2743;
        alpha = 0.1; %(mm/s)
        
        HCforce(i,1) = uf{i}(1);
        HCforce(i,2) = uf{i}(2);
    end

end