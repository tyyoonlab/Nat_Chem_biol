%% Matlab package for analysis in Magnetic Tweezers (MT) measurements %%
% % % ------------------------------------------------------------------------------------
% % %  Hyun-Kyu Choi, Minju Shon in Tae-Young Yoon's lab, SNU (Seoul National Univ)
% % % ====================================================================================
% % % ----------------------------------------------
% % % Introdcution of this script
% % % ----------------------------------------------
% % % Basically, file formats are mainly four parts: 
% % % [1] .fps (acquisition rate) [2] .cal (for off-center) [3] .r--- (Bead's position in pixel level [4] .s--- Magnet position)
% % % These componenets are originally from Labview program (by MJ) during expermient
% Note, you might have three kinds experiemnts:: Force-jump,ramp,const
%       you might be asked choose analyzing type.
% % % This script is encoded in "chronological order for analyzing"
% % % So there is no local-section/function, but local-script for each analysis
%% Sub-Script list are as below
% closest.m/correctFEC.m  = For correction of bead's off-center)
% analysis_Rm1.m = Collect whole data to single one array: x,y,z,F respectively
%  -- Note, analysis_Rm1 file has to be located in same folder with input-raw data
% eWLC_inv.m/extension_eWLC_xxx.m = For calculation of expected extension
% HSMT_forcejump.m = related analysis such as cutting traces to analze HMM ..
% HSMT_forceramp.m = crelated analysis such as FEC analysis 
% Extdistb.m = for histogram analysis at constant force rightafter/before force jump
% Modulforce.m = Linear interpretation of force-raw data due to distrupt force data points as results of inverse function
% So on....
              %%      Enjoy your analysis    %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %------------------------------------%
                %                                    %
                %              |/ / / /              %
                %          //  -    -   //           %
                %             (  @ @  )              %
                % +---------oOOo- & -oOOo--------+   %
                % |          Good Luck           |   %
                % |            H-K.C             |   %
                % +-------------------Oooo------+    %
                %           oooO      (    )         %
                %          (    )      )  /          %
                %            )  (     (__/           %
                %            (__/                    %
                %------------------------------------%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data reset
Dataclearing = 1; % Data clearing before analysis
if Dataclearing
    clear all;
    close all;
    clc;
end
%% Set temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Celcius = 23; % Degree of Celcius, ¨¬C  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find data file and collect whole data to single one array: x,y,z,F respectively

cpath = 'D:\BACKUP THIS FOLDER-190106\Experiments-200229\MT\membrane protein folding\200123 same as 2020.1.3';
dpath = 'D:\BACKUP THIS FOLDER-190106\Experiments-200229\MT\Beta2-AR\Beta2-AR\2020.2.19 same as 200208\S1\Bead1';
cd(dpath);

Dataloading = 1; % Data re-loading
if Dataloading
    display(['Loading data']);
    analysis_Rm1;
end

cd(cpath);
%% Off-center correction

Fcal = 5; % Force value during cal.
if Roff{1}{1} > 1400
    Roff{1}{1} = 400; % Our measurement bias
else
    Roff{1}{1} = Roff{1}{1};
end
disp(['%%%% Correction of Bead off-center %%%%'])
for i = 1:length(Tracez)
    Tracez_cor{i} = correctFEC(Force{i},Tracez{i},Tracex{i},Fcal,Rbead,Roff{i}{1},ori{i}{1});
end

%% Expected eWLC model calculation

Fdat_model = (0:.05:50)';

%%% Quatifying structure in polymer physics aspect %%%
extension_eWLC_b2AR; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correct quantities as adjusting your pulling assay
% Here, Ref is two 512 dsDNA handles, some polypeptide linkers and 5kMW PEG

display(['Calculate expected extension']);
T = Temp;
% Linker btw MP and DNA
Lp_PP = 0.8; Lo_AA = 0.36;
N_linker_len = 10; %10 for b2AR, 8 for GLUT3 construct
C_linker_len = 25; %25 for b2AR, 28 for GLUT3
Lo_PP = (N_linker_len+C_linker_len)*Lo_AA;
Ko_PP = 50*1000000; % polypeptide elastic modulus
zdat_model_PP = eWLC_inv(Fdat_model,Lo_PP,Lp_PP,Ko_PP,T,1);

% PEG
Lo_PEG = 570*(5/80); Lp_PEG = 0.47; % bPEG (5 kD)
Ko_PP = 50*1000000;
zdat_model_PEG = eWLC_inv(Fdat_model,Lo_PEG,Lp_PEG,Ko_PP,T,1);

% Two DNA handles
Lp_dsDNA =40.5; Ko_dsDNA = 300;
Lo_dsDNA = (2*510)*.338; % dsDNA
zdat_model_dsDNA = eWLC_inv(Fdat_model,Lo_dsDNA,Lp_dsDNA,Ko_dsDNA,T,1);

% Total extension
zdat_model_total = zdat_model_PEG + zdat_model_PP + zdat_model_dsDNA;
zdat_model_total = zdat_model_total - zdat_model_total(find(Fdat_model==Fcal));

zppGLUT3 = xp(:,1);

zdat_model_Nstate =  zdat_model_total;
zdat_model_plusMP_pp = zdat_model_total + zppGLUT3 - d_N;
zdat_model_plusMP_hs = zdat_model_total + hsxGLUT3 - d_N;

%% Intermediate states (Roughly cal)
Iu_ratio = [0.0518 0.1612 0.2705 0.3813 0.4561 0.5871 0.6561 0.7324 0.8748]; %unfolding intermediates
Ir_ratio = [0.8695 0.7119 0.4348 0.1304];  %refolding intermediates
zdat_model_plusMP_interU = {};

for ij = 1:numel(Iu_ratio)
    zdat_model_plusMP_interU{ij} = zdat_model_Nstate + (zdat_model_plusMP_pp-zdat_model_Nstate)*Iu_ratio(ij);
end

for ij = 1:numel(Ir_ratio)
    zdat_model_plusMP_interR{ij} = zdat_model_Nstate + (zdat_model_plusMP_hs-zdat_model_Nstate)*Ir_ratio(ij);
end

%% Show raw trace %%
disp(['%%%%%%%%%%%%%%%%%%% ',num2str(mean(fps{1})),'-Hz Recording traces %%%%%%%%%%%%%%%%%%%%%%'])
%Median filter set%
mfn = floor(mean(fps{1})./5);
Time{1} = (1:length(Trace{1}(:,3)))./mean(fps{1});

figure(99);
figraw = figure(99);
Screensize = get(groot, 'Screensize');
set(figraw, 'Position', [300, 150, Screensize(3)*0.8, Screensize(4)*0.75 ]);
figp = get(gcf,'Position');
set(0,'DefaultFigurePosition', figp);
                
Tsmooth = medfilt1(Trace{1}(:,3),mfn);
[axeshandles,rawforce,rawtrace]=plotyy(Time{1},Force{1},Time{1},Trace{1}(:,3));
grid on;  hold on;
set(rawforce,'Color','r','LineWidth',0.6);
set(rawtrace,'Color','k','LineWidth',1.2);

axes(axeshandles(1));
set(gca,'ytick',0:4:50);
ylabel('Force(pN)','FontSize',16);
axis([0 length(Time{1})/mean(fps{1})+1 0 50]) %max(Force{1})+4
axes(axeshandles(2));
hold on;
plot(Time{1},Tsmooth,'y','LineWidth',0.6)
ylabel('Full trace(nm)','FontSize',16)
set(gca,'ytick',-160:20:260,'linewidth', 2,'TickLength',[0.05 0.05]);
axis([0 length(Time{1})/mean(fps{1})+1 -160 260])

%% Set analysis type %%
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%  Choose Analysis  %%%%%%%%%%%%%%%%%%%%%%%%']) 
disp(['[Click] <-- : Sectioning in force-jump/constant cycle ?'])
disp(['[Click] --> : Sectioning in force-ramp cycle ?'])
disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']) 
stateanalysis_seq1 = getkey;
%% Start analysis! %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force-spectroscopy experiments are esentially divided as the two types
% Non-eq measurement : Force-jump/ Force-ramp
% Eq measurrment : Force-constant
% [ <-- ] -- set Force jump/constant to get kinetics for mapping energylandscape, extension histogram for state assignment and so on
% [ --> ] -- set Force ramp to get Force-Extension Curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if stateanalysis_seq1 == 28
    HSMT_forcejump;
elseif stateanalysis_seq1 == 29
    HSMT_forceramp;
else
    return
end
disp("complete");
%% saving files for HMM
if 1
    RawLRT = {};
    RawHRT = {};
    epath =  cpath;
    cd(epath);
    sss = find(dpath=='\');
    reforce=round(reforce.*10)/10
    Lowr = find(reforce == 4)
    Highr = find(reforce == 5)

    RawLRT{i}= Refoldingstep{Lowr}
    temp = RawLRT{i};
    save( [dpath(sss(6)+1:sss(6)+10),'_4pN_',dpath(sss(end)+1:end),' _trace',num2str(i),'.txt'],'temp','-ascii')

    
    for i = 1:length(Highr)
        RawHRT{i}= Refoldingstep{Highr(i)}
        temp = RawHRT{i};
        save( [dpath(sss(6)+1:sss(6)+10),'_5pN_',dpath(sss(end)+1:end),' _trace',num2str(i),'.txt'],'temp','-ascii')
    end
    cd(cpath);
end
%% Saving files for Hist in FECs
if 0
    fpath = cpath;
    cd(fpath);
    sss = find(dpath=='\');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_FEtrace','.mat'],'FEtracez');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_FEtrace_cor','.mat'],'FEtracez_cor');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_FETracex','.mat'],'FETracex');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_FETracey','.mat'],'FETracey');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_FEforce','.mat'],'FEforce');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_FEforce_cor','.mat'],'FEforcecorr');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_Roff','.mat'],'Roff');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_fps','.mat'],'fps');
    save([dpath(sss(6)+1:sss(6)+10),'-',dpath(sss(end)+1:end),'_ori','.mat'],'ori');
end
