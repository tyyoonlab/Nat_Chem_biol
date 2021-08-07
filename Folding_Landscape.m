%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Hyun-Kyu Choi, Hyun Gyu Kim in Tae-Young Yoon's lab % %
%------------------------------------%
%                                    %
%               ||||||               %
%              ( o  o )              %
%             (   ..   )             %
% +---------oOOo- & -oOOo--------+   %
% |        Dive SMFS analysis    |   %
% +--------------------Oooo------+   %
%           oooO     /(    )         %
%          (    )      )  /          %
%            )  (     (__/           %
%            (__/                    %
%------------------------------------%
%% MATLAB script for calculating energy landscape from constant-force ensembles
% % % This code calculates energy landscape from constant-force ensemble experiments
% % % Force-jump measurements can be analyzed by HSMT_HK.m(HSMT_forcejump.m)
%%%%%%%%%%%%%%%%%%%%% Overall code workflow %%%%%%%%%%%%%%%%%%%%%%%
% 1. Trim the constant force range in force-jump experiment
% (Can be done by HSMT_HK.m & Off-center correction can be done)
% 2. Save Trimmed z-trace as text file(off-center corrected)
% 3. Use that text file as an input
% 4. Getting z-distribution
% 5. Fourier transform and divide it by Probability of bead and handles
% 6. Calculate the probability(by sample, fourier-transformed)
% 7. Get the Free energy landscape G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
% 1. Hinczewski et al.(2013), PNAS, 110(12), 4500-4505.
% 2. Seol et al. (2007), Biophysics J., 93(12), 4360-4373.
% 3. Shon et al.(2019), Sci. Adv., 5(6), eaav1697
%% Sub-scripts lists are below
% % Main Sub-scripts
% Zh_fwlc.m ; Getting partition function using FWLC model
% FWLC model - Bending,sretching and extensiblility of biopolymer
% (script modified from Minju Shon's code)
% Zh_fjc.m ; Getting partition function using FJC model
% % Sub-scripts for Z_h_fwlc.m
% lgwt.m : Getting Legendre-Gauss nodes and weights (by Greg von Winckel)
%% Data reset
Dataclearing = 1; % Data clearing before analysis
if Dataclearing == 1
    clear all;
    close all;
    clc;
end
%%  Analysis Options
FilterResult = 1; % Filter z-trace to desirable frequency
Correct_Offcenter = 0; % Correct off-center effects
Correct_Camera = 0; % Correct camera effects
SaveResult = 0; % Save results
%% Set input-Parameters for analysis %%%
% Temperature parameters %
Celcius = 23.5; %temperature(in celcius)
Temp = 273.15 + Celcius; % Absolute temperature, K
kbT = 4.114*Temp./298.15 ; % pN*nm
beta_ = 1./kbT; % Thermodynamic beta (to distinguish from matlab beta function)

% Experimental parameters for analysis %
F0 = 5; % Force used in force-jump experiments; pN
fps = 1200; % Original frequency; Hz
Rb = 1400; % Bead radius; nm
kmax = 5; % Maximum number of data for machine learning;

Precondition = str2num(input('Filtering process ? please, input median window size (Hz) [raw = 1200Hz] : ','s'));
disp(['Cutting the trajectory? [Yes] = press Enter | [No] = press ESC ','s'])
Precondition0 = getkey;

if Correct_Offcenter == 1
    kz = 0.48; % pN/nm ; Brownian fluctuation of the magnetic bead when tethered protein is in N-state
    lp_dsDNA = 40.5; % Persistence length of DNA
    K_dsDNA = 520; % Stretch modulus of DNA
else
    kz = 0.27; % pN/nm ; Brownian fluctuation of the magnetic bead when tethered protein is in N-state
    lp_dsDNA = 11; % Persistence length of DNA
    K_dsDNA = 520; % Stretch modulus of DNA
end
if FilterResult == 1
    deltaz = 2; % 2.1 nm In 5-Hz window, Brownian fluctuation of the magnetic bead when tethered protein is in N-state
end
B = 0.25; % Normalization factor regarding x,y in probability density
exp_prefactor = exp(-0.5*beta_*kz*deltaz^2); % exp(-0.5*beta*k_z*<deltaz^2>)
% Macroscopic parameters of biopolymer like Worm-like chain model%
% PEG (biotin-PEG, 5kDa)
Lc_PEG = 570*(5/80); % Contour length of 5kDa biotin-PEG; nm
lp_PEG = 0.47; % Persistence length of 5kDa biotin-PEG; nm
K_PEG = 50*1000000; % Stretch modulus of 5kDa biotin-PEG; pN
% peptide linker
PP_len_1 = 28; % polypeptide linker 1 length(in amino acids)
PP_len_2 = 28; % polypeptide linker 2 length(in amino acids)
Lo_PP = 0.36; % Contour length of one amino-acid; nm
Lc_PP_1 = PP_len_1 * Lo_PP; % Contour length of polypeptide linker 1; nm
Lc_PP_2 = PP_len_2 * Lo_PP; % Contour length of polypeptide linker 2; nm
lp_PP = 0.39; % Persistance length of one amino-acid; nm
K_PP = 50*1000000;% Stretch modulus of amino-acid; pN
% Double-stranded DNA
dsDNA_len_1 = 512; % DNA handle (attached to bead) length; bp
dsDNA_len_2 = 512; % DNA handle (attached to surface) length; bp
Lo_dsDNA = 0.338; % Contour length of one dsDNA base-pair; nm
Lc_dsDNA_1 = dsDNA_len_1*Lo_dsDNA; % contour length of handle 1
Lc_dsDNA_2 = dsDNA_len_2*Lo_dsDNA; % contour length of handle 2
% Single-stranded DNA(used only in DNA hairpin measurements)
ssDNA_len_1 = 0; % ssDNA handle 1 length(close to bead); nt
ssDNA_len_2 = 0; % ssDNA handle 2 length(closer to surface); nt
lp_ssDNA = 0.75; % Persistence length(=kuhn length/2) of ssDNA; nm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading Data file
dpath = 'C:\Users\Owner\Desktop\Glut3_CFTR\Data and analysis\Free-energy\Test set\GlpG';
cpath = 'C:\Users\Owner\Desktop\Glut3_CFTR\Data and analysis\Free-energy';
cd(dpath);
rawdata = dir('*.txt');
S = struct(rawdata);
numdata= size(S,1);
dataraw={};

for i=1:numdata
    dataraw{i} = load(S(i).name)';
end

for i=1:numdata
    ztrace{i} = dataraw{i};
end

cd(cpath);

%% Theoretical probability density for bead, handles components
%% Calculating Pb_FT(For Bead position)
f = F0*beta_;
% Probability distribution in space domain
ksb =[-5:0.01:5];
Lb = numel(ksb); % Signal length

%Pb = exp(f.*zb).*dirac(zb-Rb)./(2*sinh(f*Rb*pi/180));
%idx = Pb  == Inf; % find Inf
%Pb(idx) = 1;     % set Inf to finite value

Pb_FT = f*sinh((f-1j*ksb).*Rb*pi/180)./((f-1j*ksb).*sinh(f*Rb*pi/180));

figure(100); hold on;
plot(ksb,Pb_FT,'k');
plot(ksb,abs((Pb_FT)),'b');
xlabel('Wave vector(k)')
ylabel('|P(k)|')
%% Calculating Ph_FT(For handle components)
% Zh_fwlc(f,Lc,lp, K,beta_,Rb,constr,fluct)
% Ph_dh1, dh2: dsDNA handles (Bead, Surface)
% Ph_peg: PEG
% Ph_ph1, ph2 : peptide handles (Bead, Surface)
% Ph_sh1, sh2 : ssDNA linkers (Bead, Surface)
% Initializing

f_rescale = f;
ks_rescale = ksb;

Ph_dh1_FT = ones(1,numel(ks_rescale));
Ph_dh2_FT = ones(1,numel(ks_rescale));
Ph_peg_FT = ones(1,numel(ks_rescale));
Ph_ph1_FT = ones(1,numel(ks_rescale));
Ph_ph2_FT = ones(1,numel(ks_rescale));
Ph_sh1_FT = ones(1,numel(ks_rescale));
Ph_sh2_FT = ones(1,numel(ks_rescale));

% PEG: half-constrained/unconstrained, no bead fluctuation
Ph_peg_FT = Zh_fwlc((f_rescale-1j*ks_rescale), Lc_PEG, lp_PEG, K_PEG, beta_, Rb, [0 1], [0 0])./Zh_fwlc(f_rescale, Lc_PEG, lp_PEG, K_PEG, beta_, Rb, [0 1], [0 0]);
% dsDNA handle 1: half-constrained/unconstrained, bead fluctuation at half-constrained end
Ph_dh1_FT = Zh_fwlc((f_rescale-1j*ks_rescale), Lc_dsDNA_1, lp_dsDNA, K_dsDNA, beta_, Rb, [1 0], [0 0])./Zh_fwlc(f_rescale, Lc_dsDNA_1, lp_dsDNA, K_dsDNA, beta_, Rb, [1 0], [0 0]);


% dsDNA handle 2
if Lc_dsDNA_2 ~= 0
    Ph_dh2_FT = Zh_fwlc((f_rescale-1j*ks_rescale), Lc_dsDNA_2, lp_dsDNA, K_dsDNA, beta_, Rb, [0 0], [0 0])./Zh_fwlc(f_rescale, Lc_dsDNA_2, lp_dsDNA, K_dsDNA, beta_, Rb, [0 0], [0 0]);
end

% % Other handle component except dsDNA_B: all unconstrained & no bead effect
% Peptide handle 1
if Lc_PP_1 ~= 0
    Ph_ph1_FT = Zh_fwlc((f_rescale-1j*ks_rescale), Lc_PP_1, lp_PP, K_PP, beta_, Rb, [0 0], [0 0])./Zh_fwlc(f_rescale, Lc_PP_1, lp_PP, K_PP, beta_, Rb, [0 0], [0 0]);
end

% Peptide handle 2
if Lc_PP_2 ~= 0
    Ph_ph2_FT = Zh_fwlc((f_rescale-1j*ks_rescale), Lc_PP_2, lp_PP, K_PP, beta_, Rb, [0 0], [0 0])./Zh_fwlc(f_rescale, Lc_PP_2, lp_PP, K_PP, beta_, Rb, [0 0], [0 0]);
end

% ssDNA handle 1
if ssDNA_len_1 ~=0
    Ph_sh1_FT = Zh_fjc((f_rescale-1j*ks_rescale), ssDNA_len_1, lp_ssDNA)./Zh_fjc(f_rescale, ssDNA_len_1, lp_ssDNA);
end

% ssDNA handle 2
if ssDNA_len_2 ~=0
    Ph_sh2_FT = Zh_fjc((f_rescale-1j*ks_rescale), ssDNA_len_2, lp_ssDNA)./Zh_fjc(f_rescale, ssDNA_len_2, lp_ssDNA);
end

% Total handle
Ph_FT = Ph_dh1_FT.*Ph_dh2_FT.*Ph_peg_FT.*Ph_ph1_FT.*Ph_ph2_FT.*Ph_sh1_FT.*Ph_sh2_FT;

figure(200);
subplot(3,3,1); hold on;
plot(ks_rescale,abs(Ph_peg_FT),'c');
plot(ks_rescale,Ph_peg_FT,'k');
subplot(3,3,2); hold on;
plot(ks_rescale,abs(Ph_dh1_FT),'color',[0 1 0]);
plot(ks_rescale,Ph_dh1_FT,'k');
subplot(3,3,3); hold on;
plot(ks_rescale,abs(Ph_dh2_FT),'color',[0.3 1 0.3]);
plot(ks_rescale,Ph_dh2_FT,'k');
subplot(3,3,4); hold on;
plot(ks_rescale,abs(Ph_ph1_FT),'r');
plot(ks_rescale,Ph_ph1_FT,'k');
subplot(3,3,5); hold on;
plot(ks_rescale,abs(Ph_ph2_FT),'color',[1 0.3 0.3]);
plot(ks_rescale,Ph_ph2_FT,'k');
subplot(3,3,7); hold on;
plot(ks_rescale,abs(Ph_sh1_FT),'b');
plot(ks_rescale,Ph_sh1_FT,'k');
subplot(3,3,8); hold on;
plot(ks_rescale,abs(Ph_sh2_FT),'color',[0.3 0.3 1]);
plot(ks_rescale,Ph_sh2_FT,'k');
subplot(3,3,6:3:9); hold on;
plot(ks_rescale,Ph_FT,'k','LineWidth',2);
plot(ks_rescale,abs(Ph_FT),'k--','LineWidth',2);

% Total effects
Phb_FT = Ph_FT.*Pb_FT;

% Space distributions by inverse Fourier Transformation
Ph = abs((ifft((Ph_FT))));
Pb = abs(ifftshift(ifft((Pb_FT))));
Phb = abs((ifft((Phb_FT))));


figure(201);
subplot(1,3,1); hold on;
plot(ks_rescale,Ph_FT,'k','LineWidth',2);
plot(ks_rescale,abs(Ph_FT),'k--','LineWidth',2);
subplot(1,3,2); hold on;
plot(ks_rescale,Pb_FT,'b','LineWidth',2);
plot(ks_rescale,abs(Pb_FT),'b--','LineWidth',2);
subplot(1,3,3); hold on;
plot(ks_rescale,Phb_FT,'g','LineWidth',2);
plot(ks_rescale,abs(Phb_FT),'g--','LineWidth',2);

weight_p = {};
mean_p = {};
sig_p = {};
zfit_p = {};
for i = 1:numdata
    %% Calculating P(z;F) & fourier transform
    % Calculating P from z-trace
    if Correct_Camera == 1
        fNyq   =  fps/2;
        % IMPLEMENTATION NEEDED
        [f_t,FTraw,T] = calc_powersp(ztrace{i},fps);
        Nf = length(f_t);
        Cblur = sinc(f_t/fps);
        FTcorr = FTraw./ Cblur;
        ind = find(f_t >= fNyq);
        FTcorr(ind) = FTraw(ind);
        
        ztrace{i} = flip(abs(ifft((FTcorr*fps))));
        
        if Precondition0 == 13
            figure(1); hold on;
            plot((1:length(ztrace{i}))./fps,ztrace{i},'k');
            plot((1:length(ztrace{i}))./fps,medfilt1(ztrace{i},round(fps/Precondition)),'r');
            disp(['Cutting trajectory for removing unwanted region']);
            pause;
            ans = getpts;
            tis0 = ans;
            ztrace_cut{i} = ztrace{i}(round(tis0(1)*fps):round(tis0(end)*fps));
            close(1);
        else
            ztrace_cut{i} = ztrace{i};
        end
        
        if FilterResult == 1
            ztrace_corr_mf{i} = medfilt1(ztrace_cut{i},round(fps/Precondition));
            ztrace_corr{i} =  ztrace_corr_mf{i}(round(fps/(2*Precondition)):end-round(fps/(2*Precondition)));
        else
            ztrace_corr{i} = ztrace_cut{i}(round(fps/(2*Precondition)):end-round(fps/(2*Precondition)));
        end
        
        figure(1); hold on;
        plot((1:length(ztrace_corr{i}))./fps,ztrace_corr{i},'b');  
        disp(['Information of z-vlaue of N-state']);
        
        pause;
        ans = getpts;
        tis = ans;
        znorm_N = mean(ztrace_corr_mf{i}(round(tis(1)*fps):round(tis(end)*fps)));
        zinitial = mean(ztrace_corr{i});
        
    else
        if Precondition0 == 13
            figure(1); hold on;
            plot((1:length(ztrace{i}))./fps,ztrace{i},'k');
            plot((1:length(ztrace{i}))./fps,medfilt1(ztrace{i},round(fps/Precondition)),'r');
            disp(['Cutting trajectory for removing unwanted region']);
            pause;
            ans = getpts;
            tis0 = ans;
            ztrace_cut{i} = ztrace{i}(round(tis0(1)*fps):round(tis0(end)*fps));
            close(1);
        else
            ztrace_cut{i} = ztrace{i};
        end
        
        if FilterResult == 1
            ztrace_corr_mf{i} = medfilt1(ztrace_cut{i},round(fps/Precondition));
            ztrace_corr{i} =  ztrace_corr_mf{i}(round(fps/(2*Precondition)):end-round(fps/(2*Precondition)));
        else
            ztrace_corr{i} = ztrace_cut{i}(round(fps/(2*Precondition)):end-round(fps/(2*Precondition)));
        end
        
        figure(1); hold on;
        plot((1:length(ztrace_corr{i}))./fps,ztrace_corr{i},'b'); 
        disp(['Information of z-vlaue of N-state']);
        pause;
        ans = getpts;
        tis = ans;
        znorm_N = mean(ztrace_corr_mf{i}(round(tis(1)*fps):round(tis(end)*fps)));
        zinitial = mean(ztrace_corr{i});
    end
    
    % Total probability density
    figure(2)
    subplot(2,2,1)
    rng default
    h = histogram(ztrace_corr{i},'Normalization','probability');
    subplot(2,2,3)
    if h.BinWidth < 0.3 || h.BinWidth > 0.5 
        numbin = round((h.BinEdges(end)-h.BinEdges(1))./0.333);
        if mod(numbin,2) == 1
            numbin = numbin - 1;
        else
            numbin = numbin;
        end
        h = histogram(ztrace_corr{i},numbin,'Normalization','probability');
        probz{i} = h.Values;
        for ii = 1:length(probz{i})
            zraw0{i}(ii) = h.BinEdges(ii)+h.BinWidth./2;
        end
    else
        numbin= h.NumBins;
        if mod(numbin,2) == 1
            numbin = numbin - 1;
        else
            numbin = numbin;
        end
        h = histogram(ztrace_corr{i},numbin,'Normalization','probability');
        z_fps = h.BinWidth;
        probz{i} = h.Values;
        for ii = 1:length(probz{i})
            zraw0{i}(ii) = h.BinEdges(ii)+h.BinWidth./2;
        end
    end
    z_fps = mean(diff(zraw0{i}));
    
    zraw{i} = zraw0{i} - zinitial; % Rearragnement of z of N state;
    zraw_N{i} = zraw0{i} - znorm_N; % Rearragnement to mean z;
    z_rearragne_ini(i) = zinitial;
    z_rearragne_N(i) = znorm_N;
    
    indn0 = probz{i} < 10^(-6) ;
    reprobz = sort(probz{i});
    indre0 = find(reprobz ~= 0);
    probz{i}(indn0) = reprobz(indre0(1)); % Inhibition for divergence
    
    subplot(2,2,2:2:4)
    plot(zraw{i},probz{i},'r');
    figure(2);
    windowWidth = 2;
    kernel = ones(windowWidth,1)/windowWidth;
    probz_lowp{i} = filter(kernel,1,probz{i});
    subplot(2,2,2:2:4); hold on;
    plot(zraw{i},probz_lowp{i},'g');
    
    % Probability distribution in space domain
    % zraw{i} : Space vector : Domain of definition
    Pcorrect{i} = probz_lowp{i}; % Praw(a function of z)
    % If consider force effects
    
    for jj = 1: numel(zraw{i})
        Ptot{i}(jj) = B*exp_prefactor*Pcorrect{i}(jj); %*exp(beta_*F0*zraw{i}(jj));
    end
    Ptot{i}=Ptot{i}./trapz(Ptot{i});
    % Sampling space-frequency = 1./z_fps
    L = length(zraw{i}); % Signal length
    weightnumdata(i) = numel(ztrace_corr{i});
    ks = linspace(-1,1,L)./(z_fps);
    
    figure(3)
    subplot(2,1,1); hold on;
    plot(zraw{i},Pcorrect{i},'k','LineWidth',1);
    plot(zraw{i},Ptot{i},'r','LineWidth',1);
    xlim([min(zraw{i})-2 max(zraw{i})+2])
    ylim([0 0.1])
    xlabel('Extension (nm)')
    ylabel('Total probability density)')
    
    % Probability distribution in wave vector (fourier transformed)
    Ptot_bare_FT{i} = fft(Pcorrect{i});
    Ptot_FT{i} = fft(Ptot{i});
    subplot(2,1,2); hold on;
    plot(ks,fftshift(abs(Ptot_bare_FT{i})),'k')
    plot(ks,fftshift(abs(Ptot_FT{i})),'r')
    xlabel('Wave vector(k)')
    ylabel('|P(k)|')
    
    
    %% rearrange using theoretical values of both handle and bead components
    for ire = 1:numel(ks)
        ind_re(ire) = find(round(ksb*100)./100 == round(ks(ire)*100)./100);
        Phb_proj{i}(ire) =  Phb(ind_re(ire));
        Ph_proj{i}(ire) =  Ph(ind_re(ire));
        Pb_proj{i}(ire) =  Pb(ind_re(ire));
    end
    
    Ph_proj{i} = Ph_proj{i}./trapz(Ph_proj{i});
    Pb_proj{i} = Pb_proj{i}./trapz(Pb_proj{i});
    Phb_proj{i} = Phb_proj{i}./trapz(Phb_proj{i});
    
    z_proj = 1/(mean(diff(ksb))*numel(ksb))*min(diff(ind_re));
    
    % Tranlsation for maximum peak to z=0;
    ind0_hb = find(Phb_proj{i} == max(Phb_proj{i}));
    ind0 = find(round(zraw{i}) == 0);
    ind0 = ind0(round(numel(ind0)/2));
    ind_shift = ind0 - ind0_hb;
    if ind_shift > 0
        Phb_proj_trans{i}(1:ind_shift) = Phb_proj{i}(L-ind_shift+1:L);
        Phb_proj_trans{i}(ind_shift+1:L) = Phb_proj{i}(1:L-ind_shift);
    else
        Phb_proj_trans{i}(1:ind0) = Phb_proj{i}(abs(ind_shift)+1:ind0_hb);
        Phb_proj_trans{i}(ind0+1:L-abs(ind_shift)) = Phb_proj{i}(ind0_hb+1:L);
        Phb_proj_trans{i}(L-abs(ind_shift)+1:L) = Phb_proj{i}(abs(ind_shift):-1:1);
    end
 
    zb{i} = linspace(Rb- z_fps*(ind0_hb-1),Rb+ z_fps*(L-ind0_hb),L);
    %zh{i} = linspace(-z_fps*(ind0_hb-1), z_fps*(L-ind0_hb),L);
    
    figure(4)
    subplot(1,3,1); hold on;
    plot(zraw{i},Ph_proj{i},'k','LineWidth',2);
    subplot(1,3,2); hold on;
    plot(zb{i},Pb_proj{i},'b','LineWidth',2);
    subplot(1,3,3); hold on;
    plot(zraw{i},Phb_proj{i},'g','LineWidth',2);
    plot(zraw{i},Phb_proj_trans{i},'r','LineWidth',2);

    %% Deconvolution in space domain
    % Gaussian fitting for both total and bead-handle systems
    
    % Total system
    options = statset('MaxIter',500);
    
    
    for k = 1:kmax
        fit_t_trial{i} =fitgmdist(ztrace_corr{i}'-zinitial,k,'Options',options,'CovarianceType','diagonal','RegularizationValue',0.1);
        AIC{i}(k)= fit_t_trial{i}.AIC;
    end
    [minAIC,numComponents] = min(AIC{i});
    AIC{i};
    Precondition1(i) = numComponents
    
    fit_tot{i} = fitgmdist(ztrace_corr{i}'-zinitial,numComponents,'Options',options,'CovarianceType','diagonal','RegularizationValue',0.1);
    p = fit_tot{i}.ComponentProportion;
    zgrid_tot{i} = zraw{i};
    for j = 1:numComponents
        nfit_tot{i}{j} = makedist('Normal','mu',fit_tot{i}.mu(j),'sigma',sqrt(fit_tot{i}.Sigma(j)));
        zfit_tot{i}{j}  = pdf(nfit_tot{i}{j},zgrid_tot{i});
        
        zfit_tot{i}{j}=  p(j).*zfit_tot{i}{j}./trapz(zfit_tot{i}{j});
        
        weight_tot{i}(j) = p(j);
        mean_tot{i}(j) = fit_tot{i}.mu(j);
        sig_tot{i}(j) = sqrt(fit_tot{i}.Sigma(j));
    end
    
    % Bead-handle system
    N_hb =1;
    zgrid_hb{i} = zraw{i};
    fit_hb = fit(zgrid_hb{i}',Phb_proj_trans{i}','gauss1');
    nfit_hb_res{i} = coeffvalues(fit_hb);
    divcheck = 0;
    for j = 1:N_hb
        if nfit_hb_res{i}(3*j-2) > 0
            zfit_hb{i}{j} = nfit_hb_res{i}(3*j-2).*exp(-((zgrid_hb{i}-nfit_hb_res{i}(3*j-1))./nfit_hb_res{i}(3*j)).^2);
            mean_hb{i}(j) = nfit_hb_res{i}(3*j-1);
            sig_hb{i}(j) = nfit_hb_res{i}(3*j)./sqrt(2);
            weight_hb{i}(j) = nfit_hb_res{i}(3*j-2);
            divcheck = divcheck+1;
        else
            zfit_hb{i}{j} = linspace(0,0,numel(zraw{i}));
            mean_hb{i}(j) = 0;
            sig_hb{i}(j) = 0;
            nfit_hb_res{i}(3*j) = 0;
            nfit_hb_res{i}(3*j-1) = 0;
            nfit_hb_res{i}(3*j-2) = 0;
            divcheck = divcheck;
        end
    end
    
    if divcheck ~= N_hb
        nfit_hb_res = {};
        zfit_hb{i} = {};
        mean_hb{i} = 0;
        sig_hb{i} = 0;
        
        N_hb =1;
        fit_hb = fit(zgrid_hb{i}',Phb_proj_trans{i}','gauss1');
        nfit_hb_res{i} = coeffvalues(fit_hb);
        j=1;
        zfit_hb{i}{j} = nfit_hb_res{i}(1).*exp(-((zgrid_hb{i}-nfit_hb_res{i}(2))./nfit_hb_res{i}(3)).^2);
        mean_hb{i}(j) = nfit_hb_res{i}(2);
        sig_hb{i}(j) = nfit_hb_res{i}(3)/sqrt(2);
        weight_hb{i}(j) = nfit_hb_res{i}(1);
    end
    
    weight_hb{i} = weight_hb{i}./sum(weight_hb{i});
    
    figure(5)
    subplot(1,2,1)
    hold on;
    for j = 1:numComponents
        plot(zgrid_tot{i},zfit_tot{i}{j},'LineWidth',1);
    end
    plot(zraw{i},Pcorrect{i},'k','LineWidth',2);
    plot(zraw{i},Ptot{i},'r','LineWidth',2);
    ylim([0 0.1])
    xlim([min(zraw{i})-2 max(zraw{i})+2])
    xlabel('Extension (nm)')
    ylabel('Normalized probability density of total system')
    
    subplot(1,2,2)
    plot(zraw{i},Phb_proj_trans{i},'g','LineWidth',2);
    hold on;
    for j = 1:N_hb
        plot(zgrid_hb{i},zfit_hb{i}{j},'LineWidth',1);
    end
    plot(fit_hb,zgrid_hb{i},Phb_proj_trans{i});
    xlim([min(zraw{i})-2 max(zraw{i})+2])
    ylim([0 0.3])
    xlabel('Extension (nm)')
    ylabel('Normalized probability density of bead+handle system')
    
    
    %% Get pure protein's probabilty density
    %figure(6)
    for ii = 1:Precondition1(i)
        weight_p{i}(ii) = weight_tot{i}(ii);
        mean_corr_sum = 0;
        
        for jj = 1:N_hb
            mean_corr_sum = mean_corr_sum + weight_hb{i}(jj).*mean_hb{i}(jj);
        end
        
        mean_p{i}(ii) = mean_tot{i}(ii)  - mean_corr_sum;
    end
    
    for ii = 1:Precondition1(i)
        % Deconvolution of variance using covariance
        if N_hb == 1
            rnd_tothb = normrnd(mean_tot{i}(ii),sig_tot{i}(ii),[1 1000*L]).*normrnd(mean_hb{i}(N_hb),sig_hb{i}(N_hb),[1 1000*L]);
            E_tothb =  (1/(1000*L))*sum(rnd_tothb,2);
            E_hb = mean_hb{i}(N_hb);
            E_tot = mean_tot{i}(ii);     
            sig_p{i}(ii) = sig_tot{i}(ii).^2  - sig_hb{i}(N_hb).^2 - 2*(E_tothb-E_hb*E_tot);
        else      
            sig_corr_sum = 0;
            for k = 1:N_hb
                rnd_tothb = normrnd(mean_tot{i}(ii),sig_tot{i}(ii),[1 1000*L]).*normrnd(mean_hb{i}(k),sig_hb{i}(k),[1 1000*L]);
                E_tothb =  (1/(1000*L))*sum(rnd_tothb,2);
                E_hb = mean_hb{i}(k);
                E_tot = mean_tot{i}(ii);
                sig_corr_sum = sig_corr_sum + sig_hb{i}(k).^2 + 2*(E_tothb-E_hb*E_tot);
            end
            sig_p{i}(ii) = sig_tot{i}(ii).^2  - sig_corr_sum;
        end
    end
    sig_p{i} = sqrt(sig_p{i});
    %% Free energy landscape
    
    for j =  1:Precondition1(i)
        zfit_p{i}{j} = exp(-((zraw{i}-mean_p{i}(j))./sig_p{i}(j)).^2);
        zfit_p{i}{j} = zfit_p{i}{j}./trapz(zfit_p{i}{j});
    end
    
    for i_array = 1:L
        pz_mp = 0;
        for i_run = 1:Precondition1(i)
            pz_mp = pz_mp+weight_p{i}(i_run)*zfit_p{i}{i_run}(i_array);
        end
        Pmp{i}(i_array)  = pz_mp;
    end
    
    figure(6); hold on;
    plot(zraw{i},Pcorrect{i},'k','LineWidth',1);
    plot(zraw{i},Ptot{i},'r','LineWidth',1.5);
    plot(zraw{i},Pmp{i},'b','LineWidth',2);
    ylim([0 max(Pmp{i})+0.05])
    xlim([min(zraw{i})-2 max(zraw{i})+2])
    xlabel('Extension (nm)')
    ylabel('Normalized probability density of total system')
    
    G_bare{i} = -log(abs(ifft(fftshift(Ptot_bare_FT{i}))))./log(exp(1));
    G_tot{i} = -log(abs(ifft(fftshift(Ptot_FT{i}))))./log(exp(1));
    G_corr{i} = -log(Pmp{i})./log(exp(1)); % Energy landscape
    
    delGr = + beta_*F0*abs(zraw_N{i});
    figure(7);
    subplot(1,2,1);hold on;
    plot(zraw_N{i},G_bare{i},'k','LineWidth',2);
    plot(zraw_N{i},G_tot{i},'r','LineWidth',2);
    plot(zraw_N{i},G_corr{i},'b','LineWidth',2);
    subplot(1,2,2);hold on;
    plot(zraw_N{i},G_bare{i}+delGr,'k','LineWidth',2);
    plot(zraw_N{i},G_tot{i}+delGr,'r','LineWidth',2);
    plot(zraw_N{i},G_corr{i}+delGr,'b','LineWidth',2);
    
    pause;
    close all;
end

%% Weighting Coarse-Grained Free-energy landscape

% Weighting distribution along number of data
z_protein = round(min(cell2mat(ztrace_corr)))-1:0.1:round(max(cell2mat(ztrace_corr)))+1;
for  i = 1:numdata
    for j =  1:Precondition1(i)
        zfit_p{i}{j} = exp(-((z_protein-z_rearragne_ini(i)-mean_p{i}(j))./sig_p{i}(j)).^2);
        zfit_p{i}{j} = zfit_p{i}{j}./trapz(zfit_p{i}{j});
    end
    
    for i_array = 1:numel(z_protein)
        pz_mp = 0;
        for i_run = 1:Precondition1(i)
            pz_mp = pz_mp+weight_p{i}(i_run)*zfit_p{i}{i_run}(i_array);
        end
        Pmp_repmat{i}(i_array)  = pz_mp;
    end
    G_repmat{i} = -log(Pmp_repmat{i})./log(exp(1));
    z_final{i} = z_protein - round(z_rearragne_N(i)*10)./10;
    
    G_repmat_0{i} = G_repmat{i}+beta_*F0*abs(z_final{i});
    
    indgmin_final = find(G_repmat_0{i} == min(G_repmat_0{i}));
    % Rearrange along y-axis
    G_repmat_0{i}  = G_repmat_0{i} - G_repmat_0{i} (indgmin_final);
    
    % Rearrange along x-axis
    z_final_0{i} =  z_final{i}-  z_final{i}(indgmin_final);
    z_final_0{i} = round(z_final_0{i}.*10)./10;
    ind_mp_min(i) = z_final_0{i}(1);
    ind_mp_max(i) = z_final_0{i}(end);
end

cut_start = max(ind_mp_min);
cut_end = min(ind_mp_max);
for i = 1:numdata
    a = find(z_final_0{i}==cut_start);
    b = find(z_final_0{i}==cut_end);
    z_final_0_re{i} = z_final_0{i}(a:b);
    G_repmat_0_re{i} = G_repmat_0{i}(a:b);
end

z_mp_final = cut_start:0.1:cut_end;

for j_array = 1:numel(z_mp_final)
    pz_mp_remat_weight = 0;
    pz_mp_remat_mean =  0;
    for i = 1:numdata
        pz_mp_remat_weight = pz_mp_remat_weight + weightnumdata(i)*G_repmat_0_re{i}(j_array)./sum(weightnumdata));
        pz_mp_remat_mean =  pz_mp_remat_mean + G_repmat_0_re{i}(j_array)./(numel(weightnumdata));
    end
    G_mp_final_weight(j_array) = pz_mp_remat_weight;
    G_mp_final_mean(j_array) = pz_mp_remat_mean;   
end

figure(1); hold  on;
for i = 1:numdata
    plot(z_final_0_re{i},G_repmat_0_re{i},'LineWidth',1);
    pause;
end
plot(z_mp_final,G_mp_final_weight,'r','LineWidth',3);
plot(z_mp_final,G_mp_final_mean,'b','LineWidth',3);
xlim([round(cut_start) round(cut_end)]);
ylim([-2.5 35]);
%% Calculating P & inverse fourier transform
%% %%%%%%%%%%%%% Numerically unstable %%%%%%%%%%%%%%%%
if 0
    Pmp_FT{i} = Ptot_FT{i}./fftshift(Phb_FT); % Probability distribution of sample(fourier transformed)
    Pmp{i} = abs(ifft(Pmp_FT{i}));
    Pmp{i} = Pmp{i}./trapz(Pmp{i});
    
    figure(8);
    subplot(1,2,1);hold on;
    plot(zraw{i},Pcorrect{i},'k');
    plot(zraw{i},Ptot{i},'r');
    plot(zraw{i},Pmp{i},'b');
    xlabel('Relative extension (nm)');
    ylabel('P(z)');
    subplot(1,2,2);hold on;
    plot(ks,fftshift(abs(Ptot_bare_FT{i})),'k--');
    plot(ks,fftshift(abs(Ptot_FT{i})),'r--','LineWidth',2);
    plot(ks,abs(fftshift(Pmp_FT{i})),'b--');
    plot(ks,fftshift(Ptot_bare_FT{i}),'k');
    plot(ks,fftshift(Ptot_FT{i}),'r');
    plot(ks,fftshift(Pmp_FT{i}),'b');
    xlabel('Wave vector(k,1/nm)');
    ylabel('|P(k)|');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
cd(cpath);
%% Saving results
if SaveResult == 1
    cd(epath);
    % NEEDS IMPLEMENTATION
end

