%%%%%%% Sub .m file for analyis in HSMT %%%%%%
disp(['%%%%%%%%%%%%%%%%%%%  Force-jump/constant Measurement  %%%%%%%%%%%%%%%%%%%%%%%']);
pause;
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
disp(['[Click] <-- : Trace cutting for HMM'])
disp(['[Click] --> : Trace histogram for state analysis'])
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
stateanalysis_seq2 = getkey;

tis0 ={}; % check points of unfolding trace;
tis1 = {}; % check points of whole traces;
tis2 = {}; % check points of refolding trace;
tis3 = {}; % check points for z-drift correction;
tis4 = {}; % check points of refolding steps;

Retracez = {};
Retracez_cor ={};
Resmooth = {};
ReTime = {};
Refoldingstep ={};
Stretchingstep ={};
SandRstep = {};
nc = 100; % number of analysis-iteration
waitingTime = zeros(nc,1);
reforce = zeros(nc,1);
forcejump = zeros(nc,1);
populations={};
locations={};
% Temporal cutting check-point
a1=1;
c1=1;

for i =1:100
    pause;
    disp("%%Press Any Key%%");
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    disp(['[Click] <-- : Analyzing again previous cycle'])
    disp(['[Click] --> : Analyzing next cycle'])
    disp(['[Click] esc : Stop analyzing'])
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    stateanalysis_seq0 = getkey;
    
    if stateanalysis_seq0 == 28
        i = i - 1;
        if i == 0
            i = 1;
        end
    elseif stateanalysis_seq0 == 29
        i = i;
    else
        return
    end
    
    if stateanalysis_seq2 == 28 ||  stateanalysis_seq2 == 29
        dcheck=mean(Force{1}(round(a1):round(c1)));
        figure(99);
        axeshandles = gca;
        hold on;
        plot(a1./mean(fps{1}),dcheck, 'gx', 'linewidth', 3, 'markersize', 20);
        plot(c1./mean(fps{1}),dcheck, 'gx', 'linewidth', 3, 'markersize', 20);
        
        if 1
            figure(99);
            disp("refolding section of the trace");
            disp("press enter to set the start & end. double-click to end");
            pause;
            getpts
            tis2{i} = ans;
            
            a2 = tis2{i}(1)* mean(fps{1});
            c2 = tis2{i}(end)* mean(fps{1});
            waitingTime(i) = (c2 - a2)./ mean(fps{1});
            reforce(i) = mean(Force{1}((a2):(c2)))
        end
        
        if 1
            figure(99);
            disp("unfolding section of the trace");
            disp("press enter to set the start & end. double-click to end");
            pause;
            getpts
            tis0{i} = ans;
            
            a0 = tis0{i}(1)* mean(fps{1});
            c0 = tis0{i}(end)* mean(fps{1});
            forcejump(i) = mean(Force{1}((a0):(c0)))
        end
        
        if 1
            disp("section of the trace to be analyzed");
            disp("press any key to set start & end points. click to set start point. double-click to end");
            pause;
            getpts
            tis1{i} = ans;
            
            a1 = tis1{i}(1)* mean(fps{1});
            c1 = tis1{i}(end)* mean(fps{1});
            
            Retracez{i} = Trace{1}(a1:c1,3);
            Retracez_cor{i} = Tracez_cor{1}(a1:c1);
            Retracex{i} = Trace{1}(a1:c1,1);
            Retracey{i} = Trace{1}(a1:c1,2);
            Reforce{i} =  Force{1}(a1:c1);
            ReTime{i} = (1:length(Retracez{i}))/mean(fps{1});
            
            
            figure(200+i);
            plot(ReTime{i},Retracez{i},'Color',[0.25 0.25 0.25],'LineWidth',1);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('z (nm)');
            title(['Zoom-in trace versus Time']);
        end
        
        %% Matching with theoretical extension at each states
        if 1
            % Theoretical Expected extension of MP
            
            fl = find(round(Fdat_model*10)/10==round(reforce(i).*10)/10);
            fjl = find(round(Fdat_model*10)/10==round(forcejump(i).*10)/10);
            
            z_Uz = zdat_model_plusMP_hs(fl(round(length(fl)/2)));
            z_N_fr = zdat_model_Nstate(fl(round(length(fl)/2)));
            %refoldingI1 = refoldingz + (z_Uz-refoldingz)*rfRatioI1;
            
            z_pp = zdat_model_plusMP_pp(fjl(round(length(fjl)/2)));
            z_N_fj = zdat_model_Nstate(fjl(round(length(fjl)/2)));
            
            for ij = 1:numel(Iu_ratio)
                z_pi(ij) = zdat_model_plusMP_interU{ij}(fjl(round(length(fjl)/2)));
            end
            
            for ij = 1:numel(Ir_ratio)
                z_hi(ij) = zdat_model_plusMP_interR{ij}(fl(round(length(fl)/2)));
            end
            
            display(['Shift traces for matching expected extension of cycle' num2str(i)])
            
            m = 2;
            x = 16;
            figure(200+i);
            f=figure(200+i);
            Screensize = get(groot, 'Screensize' );
            set(f, 'Position', [300, 125, Screensize(3)*0.8, Screensize(4)*0.75 ]);
            figp = get(gcf,'Position');
            set(0,'DefaultFigurePosition', figp);
            display(['Left arrow == shift all traces raw data -',num2str(x)]);
            display(['Right arrow == shift all traces raw data +',num2str(x)]);
            display(['Up arrow == shift all traces raw data +',num2str(m)]);
            display(['Down arrow == shift all traces raw data -',num2str(x)]);
            
            while 1
                f;
                axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
                plot(get(gca,'xlim'), [z_N_fr z_N_fr],'--r','LineWidth',2);
                hold on;
                plot(get(gca,'xlim'),[z_Uz z_Uz],'--m','LineWidth',2);
                hold on;
                plot(get(gca,'xlim'),[z_N_fj z_N_fj],'r','LineWidth',2);
                plot(get(gca,'xlim'),[z_pp z_pp],'--b','LineWidth',2);
                
                for ij = 1:numel(Iu_ratio)
                    %plot(get(gca,'xlim'),[z_pi(ij) z_pi(ij)],'--g','LineWidth',2);
                end
                for ij = 1:numel(Ir_ratio)
                    plot(get(gca,'xlim'),[z_hi(ij) z_hi(ij)],'--g','LineWidth',2);
                end
                
                plot(ReTime{i},medfilt1(Retracez{i},mean(fps{1})./5),'k','LineWidth',0.5);
                %plot(ReTime{i},medfilt1(Retracez_cor{i},mean(fps{1})./5),'g','LineWidth',0.5);
                axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
                grid on;
                set(gca,'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
                xlabel('Time (s)'); ylabel('z (nm)');
                title(['Zoom-in trace versus Time']);
                ans = getkey
                if ans == 28
                    Retracez{i} = Retracez{i} - m;
                    %Retracez_cor{i} = Retracez_cor{i} - m;
                elseif ans == 29
                    Retracez{i} = Retracez{i} + m;
                    %Retracez_cor{i} = Retracez_cor{i} + m;
                elseif ans == 30
                    Retracez{i} = Retracez{i} + x;
                    %Retracez_cor{i} = Retracez_cor{i} + x;
                elseif ans == 31
                    Retracez{i} =Retracez{i} - x;
                    %Retracez_cor{i} = Retracez_cor{i} - x;
                else
                    break
                end
                clf(f);
            end
            close(f);
        end
        %% Setting interval such as refolding region, unfolding region (I1,I2,U states etc)
        sz = 10;
        figure(300+i);
        axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
        plot(get(gca,'xlim'), [z_N_fr z_N_fr],'--r','LineWidth',2);
        hold on;
        plot(get(gca,'xlim'),[z_Uz z_Uz],'--m','LineWidth',2);
        %plot(get(gca,'xlim'),[z_N_fj z_N_fj],'r','LineWidth',2);
        plot(get(gca,'xlim'),[z_pp z_pp],'--b','LineWidth',2);
        for ij = 1:numel(Iu_ratio)
            %plot(get(gca,'xlim'),[z_pi(ij) z_pi(ij)],'--g','LineWidth',2);
        end
        
        for ij = 1:numel(Ir_ratio)
            plot(get(gca,'xlim'),[z_hi(ij) z_hi(ij)],'--g','LineWidth',2);
        end
        hold on;
        plot(ReTime{i},Retracez{i},'color',[0.7 0.7 0.7 0.2],'LineWidth',0.2);
        %plot(ReTime{i},Retracez_cor{i},'color',[1 0.6 0.6],'LineWidth',0.2);
        plot(ReTime{i},medfilt1(Retracez{i},mean(fps{1})./5),'k','LineWidth',0.6);
        %plot(ReTime{i},medfilt1(Retracez_cor{i},mean(fps{1})./5),'r','LineWidth',0.6);
        axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
        set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Time (s)'); ylabel('z (nm)');
        %set(gcf,'Renderer','painters')
    else
        return
    end
    
    % long-period range z-drift correction - assumed as linear drift
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    disp(['[Click] <-- : Skip correction of z-drift in linear-drift assumption'])
    disp(['[Click] anybutton : Apply correction of z-drift in linear-drift assumption'])
    disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
    stateanalysis_seq3 = getkey;
    
    if stateanalysis_seq3  == 28
        close(300+i);
        figure(300+i);
        for fii = 1:2
            if fii == 1
                vistrace = Retracez{i};
            elseif fii == 2
                vistrace = Retracez_cor{i};
            end
            subplot(1,2,fii)
            axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
            plot(get(gca,'xlim'), [z_N_fr z_N_fr],'--r','LineWidth',2);
            hold on;
            plot(get(gca,'xlim'),[z_Uz z_Uz],'--m','LineWidth',2);
            plot(get(gca,'xlim'),[z_N_fj z_N_fj],'r','LineWidth',2);
            plot(get(gca,'xlim'),[z_pp z_pp],'--b','LineWidth',2);
            for ij = 1:numel(Iu_ratio)
                plot(get(gca,'xlim'),[z_pi(ij) z_pi(ij)],'--g','LineWidth',2);
            end
            for ij = 1:numel(Ir_ratio)
                plot(get(gca,'xlim'),[z_hi(ij) z_hi(ij)],'--g','LineWidth',2);
            end
            hold on;
            plot(ReTime{i},vistrace,'color',[0.7 0.7 0.7],'LineWidth',0.2);
            plot(ReTime{i},medfilt1(vistrace,mean(fps{1})./5),'k','LineWidth',0.6);
            hold on;
            axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('z (nm)');
            %set(gcf,'Renderer','painters')
        end
    else
        figure(300+i);
        disp("press any key to set start&end points to set drift pos.");
        pause;
        getpts
        tis3{i} = ans;
        
        zc1 = round(tis3{i}(1)*mean(fps{1}));
        zc2 = round(tis3{i}(2)*mean(fps{1}));
        zc3 = round(tis3{i}(end-1)*mean(fps{1}));
        zc4 = round(tis3{i}(end)*mean(fps{1}));
        
        zdrift = (mean(Retracez{i}(zc3:zc4))-mean(Retracez{i}(zc1:zc2)))./(zc4-zc1) ;% nm/frame(#)
        for loc = 1:length(Retracez{i})
            Retracez{i}(loc) = Retracez{i}(loc) - zdrift*loc ;
            Retracez_cor{i}(loc) = Retracez_cor{i}(loc) - zdrift*loc ;
        end
        close(300+i);
        figure(300+i);
        for fii = 1:2
            if fii == 1
                vistrace = Retracez{i};
            elseif fii == 2
                vistrace = Retracez_cor{i};
            end
            subplot(1,2,fii)
            axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
            plot(get(gca,'xlim'), [z_N_fr z_N_fr],'--r','LineWidth',2);
            hold on;
            plot(get(gca,'xlim'),[z_Uz z_Uz],'--m','LineWidth',2);
            plot(get(gca,'xlim'),[z_N_fj z_N_fj],'r','LineWidth',2);
            plot(get(gca,'xlim'),[z_pp z_pp],'--b','LineWidth',2);
            for ij = 1:numel(Iu_ratio)
                plot(get(gca,'xlim'),[z_pi(ij) z_pi(ij)],'--g','LineWidth',2);
            end
            
            for ij = 1:numel(Ir_ratio)
                plot(get(gca,'xlim'),[z_hi(ij) z_hi(ij)],'--g','LineWidth',2);
            end
            hold on;
            plot(ReTime{i},vistrace,'color',[0.7 0.7 0.7],'LineWidth',0.2);
            plot(ReTime{i},medfilt1(vistrace,mean(fps{1})./5),'k','LineWidth',0.6);
            hold on;
            axis([0 length(Retracez{i})./mean(fps{1}) z_N_fr-40 z_pp+40]);
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('z (nm)');
            %set(gcf,'Renderer','painters')
        end
    end
    
    if stateanalysis_seq2 == 28
        figure(400+i);
        [axeshandles,rawtrace,rawforce]=plotyy((1:length(Retracez{i}))./mean(fps{1}),Retracez{i},(1:length(Reforce{i}))./mean(fps{1}),Reforce{i});
        grid on;
        set(rawtrace,'Color','k','LineWidth',0.5);
        set(rawforce,'Color','r','LineWidth',0.5);
        axes(axeshandles(1));
        hold on;
        plot((1:length(Retracez{i}))./mean(fps{1}),medfilt1(Retracez{i},mean(fps{1})/5),'y','LineWidth',1)
        plot((1:length(Retracez_cor{i}))./mean(fps{1}),medfilt1(Retracez_cor{i},mean(fps{1})/5),'g','LineWidth',1)
        ylabel('Refolding trace(nm)','FontSize',16)
        axis([0 length(Retracez{i})/mean(fps{1})+1 min(Retracez{i})-40 max(vistrace)+40])
        axes(axeshandles(2));
        ylabel('Force(pN)','FontSize',16);
        axis([0 length(Reforce{i})/mean(fps{1})+1 0 max(Reforce{i})+4])
        
        figure(400+i);
        pause;
        
        getpts
        tis4{i} = ans;
        arf = tis4{i}(1)*mean(fps{1});
        crf = tis4{i}(end)*mean(fps{1});
        
        Refoldingstep{i} = Retracez{i}(arf:crf);
        Refoldingstep_corr{i} = Retracez_cor{i}(arf:crf);
        
    elseif stateanalysis_seq2 == 29
        disp(['Extension distribution in unfolding phase'])
        % close(300+i);
        figure(400+i);
        [axeshandles,rawtrace,rawforce]=plotyy((1:length(Retracez{i}))./mean(fps{1}),Retracez{i},(1:length(Reforce{i}))./mean(fps{1}),Reforce{i});
        grid on;
        set(rawtrace,'Color','k','LineWidth',0.5);
        set(rawforce,'Color','r','LineWidth',0.5);
        axes(axeshandles(1));
        hold on;
        plot((1:length(Retracez{i}))./mean(fps{1}),medfilt1(Retracez{i},mean(fps{1})/5),'y','LineWidth',1)
        plot((1:length(Retracez_cor{i}))./mean(fps{1}),medfilt1(Retracez_cor{i},mean(fps{1})/5),'g','LineWidth',1)
        ylabel('Refolding trace(nm)','FontSize',16)
        axis([0 length(Retracez{i})/mean(fps{1})+1 min(Retracez{i})-40 max(vistrace)+40])
        axes(axeshandles(2));
        ylabel('Force(pN)','FontSize',16);
        axis([0 length(Reforce{i})/mean(fps{1})+1 0 max(Reforce{i})+4])
        
        pause;
        
        getpts
        tis4{i} = ans;
        
        brf = tis4{i}(1)*mean(fps{1});
        crf = tis4{i}(end)*mean(fps{1});
        
        Refoldingstep{i} = Retracez{i}(brf-0.025*mean(fps{1}):brf);
        Stretchingstep{i} = Retracez{i}(brf:crf+0.025*mean(fps{1}));
        SandRstep{i} = Retracez{i}(brf-0.025*mean(fps{1}):crf+0.025*mean(fps{1}));
        
        close(400+i);
        
        stretdata = Stretchingstep{i}; %
        nbins{i} = ceil((max(SandRstep{i})-min(SandRstep{i})));
        [counts,bins] = hist(stretdata,nbins{i}); %#get populations{i}  and bin locations
        
        populations{i} = counts ;
        locations{i} =bins;
        
        updist{i} = Retracez{i}(floor(crf):floor(crf+0.1*mean(fps{1})));
        upgfit=fitgmdist(updist{i},1);
        nfit = makedist('normal',upgfit.mu,sqrt(upgfit.Sigma));
        
        %p = upgfit.ComponentProportion;
        
        zgrid{i} = linspace(min(SandRstep{i}),max(SandRstep{i}),nbins{i}*20)';
        zfit = pdf(nfit,zgrid{i});
        zgfit{i} = ceil(max(hist(updist{i},floor(max(updist{i})-min(updist{i}))))/max(zfit))*zfit;
        
        threshold{i} = upgfit.mu - (z_pp-z_N_fj);
        
        upmean(i) = upgfit.mu;
        upsig(i) = upgfit.Sigma;
        
        if 1
            % find thresold line level
            netbins{i} = locations{i} - threshold{i};
            ceilposi = find(abs(ceil(netbins{i}))==1);
            floorposi =  find(abs(floor(netbins{i}))==1);
            bothposi = find(populations{i}==min(populations{i}));
            
            if  (isempty(floorposi) == 1) && (isempty(ceilposi) == 1)
                floorposi = bothposi;
                ceilposi = bothposi;
                
                for h = 1:length(locations{i})
                    cutshift{i}(h)= round(min(populations{i}));
                end
                
            elseif isempty(floorposi) == 1
                floorposi = ceilposi;
                binposi{i}(1) = ceilposi(1);
                binposi{i}(2) = floorposi(1);
                cutoffz = floor(sum(binposi{i})./2);
                
                for h = 1:length(locations{i})
                    cutshift{i}(h)= round(mean(populations{i}(1:cutoffz)));
                end
                
            elseif isempty(ceilposi) == 1
                ceilposi = floorposi;
                binposi{i}(1) = ceilposi(1);
                binposi{i}(2) = floorposi(1);
                cutoffz = floor(sum(binposi{i})./2);
                
                for h = 1:length(locations{i})
                    cutshift{i}(h)= round(mean(populations{i}(bothposi)));
                end
            else
                cutoffz = bothposi(1);
                for h = 1:length(locations{i})
                    cutshift{i}(h)= round(mean(populations{i}(1:cutoffz)));
                end
            end
        end
        
        figure(600+i);
        subplot(1,2,1)
        hold on; box on; grid on;
        axis([0 length(Stretchingstep{i})./mean(fps{1}) min(Stretchingstep{i})-10 max(Stretchingstep{i})+10])
        plot(get(gca,'xlim'), [upgfit.mu upgfit.mu],'b','LineWidth',2);
        plot(get(gca,'xlim'), [threshold{i} threshold{i}],'r','LineWidth',2);
        plot((1:length(Stretchingstep{i}))./mean(fps{1}),Stretchingstep{i},'color',[0.65 0.65 1],'LineWidth',1);
        plot((1:length(Stretchingstep{i}))./mean(fps{1}),medfilt1(Stretchingstep{i},mean(fps{1})/20),'b','LineWidth',1.5);
        set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Time (s)'); ylabel('z (nm)');
        
        subplot(1,2,2)
        axis([0 max(populations{i} )+5 min(Stretchingstep{i})-10 max(Stretchingstep{i})+10])
        barh(locations{i},populations{i} ,'FaceColor',[0.75 0.25 1],'EdgeColor',[0.75 0.25 1],'facealpha',0.75);
        hold on;
        barh(locations{i},cutshift{i},'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0.75 0.75 0.75],'facealpha',0.75);
        hold on;
        plot(zgfit{i},zgrid{i},'k','LineWidth',1.5);
        plot(get(gca,'xlim'), [upgfit.mu upgfit.mu],'b','LineWidth',2);
        plot(get(gca,'xlim'), [threshold{i} threshold{i}],'r','LineWidth',2);
        axis([0 max(populations{i} )+5 min(Stretchingstep{i})-10 max(Stretchingstep{i})+10])
        set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
        xlabel('Count (#)');
        
        %%%%%%%%%%%%%%%%%
        zup = 200; %%%%%%
        %%%%%%%%%%%%%%%%%
        
        zdelshift(i) = upmean(i) - zup;
        
        zgrid{i} = zgrid{i} - zdelshift(i);
        Stretchingstep{i} = Stretchingstep{i} - zdelshift(i);
        Refoldingstep{i} = Refoldingstep{i} - zdelshift(i);
        SandRstep{i} = SandRstep{i} - zdelshift(i);
        locations{i} = locations{i} - zdelshift(i);
        
        threshold_2(i) = zup - (z_pp-z_N_fj);
        
        Nz(i) = zdat_model_Nstate(find (Fdat_model == 5))- zdelshift(i);
        UHz(i) = z_Uz- zdelshift(i) - z_N_fr;
    end
end

if stateanalysis_seq2 == 29
    pause;
    % Calculate noise level
    for i = 1:nc
        coc(i) = mean(cutshift{i});
    end
    cocl = mean(coc);
    
    figure(800)
    pl = ceil(floor(sqrt(nc)*10)/10);
    for i = 1:nc
        subplot(4,2*pl,2*i-1)
        hold on; box on; grid on;
        axis([0 length(Stretchingstep{i})./mean(fps{1}) min(Stretchingstep{i})-10 max(Stretchingstep{i})+10])
        plot(get(gca,'xlim'), [zup zup],'b','LineWidth',1);
        plot(get(gca,'xlim'), [threshold_2(i) threshold_2(i)],'r','LineWidth',1);
        plot((1:length(Refoldingstep{i}))./mean(fps{1}), Refoldingstep{i},'color',[1 0.55 0.55],'LineWidth',1);
        plot((length(Refoldingstep{i})+1:length(Refoldingstep{i})+length(Stretchingstep{i}))./mean(fps{1}),Stretchingstep{i},'color',[0.55 0.55 1],'LineWidth',1);
        plot((1:length(SandRstep{i}))./mean(fps{1}),medfilt1(SandRstep{i},mean(fps{1})/20),'k','LineWidth',1.5);
        set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Time (s)'); ylabel('z (nm)');
        
        subplot(4,2*pl,2*i)
        axis([0 max(populations{i} )+5 min(Stretchingstep{i})-10 max(Stretchingstep{i})+10])
        hold on; box on; grid on;
        plot(get(gca,'xlim'), [zup zup],'b','LineWidth',2);
        plot(get(gca,'xlim'), [threshold_2(i) threshold_2(i)],'r','LineWidth',2);
        barh(locations{i},populations{i} ,'FaceColor',[0.75 0.25 1],'EdgeColor',[0.75 0.25 1],'facealpha',0.75);
        barh(locations{i},cutshift{i},'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0.75 0.75 0.75],'facealpha',0.75);
        plot(zgfit{i},zgrid{i},'k','LineWidth',1.5);
        axis([0 max(populations{i} )+5 min(Stretchingstep{i})-10 max(Stretchingstep{i})+10])
        set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
        xlabel('Count (#)');
        set(gca,'YTickLabel',[]);
    end
    
    figure(801);
    sc=55;
    subplot(1,2,1)
    axis([0 nc+1 0 15])
    plot(get(gca,'xlim'), [cocl cocl],'r','LineWidth',1);
    hold on;
    scatter(1:length(coc),coc,sc,'ko','filled')
    grid on;
    axis([0 nc+1 0 10])
    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
    ylabel('Noise threshold level (#)');
    hold off;
    
    subplot(1,2,2)
    hold on;
    scatter(reforce,forcejump,sc,'bo','filled')
    grid on;
    axis([0 10 10 25])
    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
    title('F_{f} vs F_{uf}');
    hold off;
    
    pause;
    
    if 1
        Extdistb
    end
end