Fmag =  Fdat_model'; % Exerted magnetic force
Temp = 273.15 + Celcius; % Absolute temperature, K
kbT = 4.114*Temp./298.15 ; % pN*nm

%%%%%%%%%%%%%%%%%%%% Inf.of target MP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 12 ; % number of helices in GLUT3
Np = 463;  % Number of amino-acid residues of GLUT3
lp = 0.36; % contour length of one amino-acid
Pp = 0.4; % persistance length of one amino-acid
bk = 0.55; % Kuhn length of amino-acids


lh = 0.16; % contour length of alpha helix per one monomer amino-acid
Ph = 9.17; % persistance length of one alpha helix
Kaa = 10^(12);% Stretch modulus of amino-acid; pN

f = Fmag./(kbT)';% 1/nm
% Number of residues for each helix
% GLUT3 is suspected to be unfolded from C-term(H6) to N-term(H1)
NH = zeros(m,1);
NH = [22 24 21 29 26 21 36 24 25 31 30 20]'; % Number of amino-acid in TMD from PDB

Nl = zeros(m,1);
Nl = [11 2 6 2 10 11 3 4 5 6 2 7]'; % Number of amino-acid in polypeptide form in Linker from PDB

NHl = zeros(m,1);
NHl= [15 8 0 0 0 48 0 0 0 0 0 14]'; % Number of amino-acid in helix form in Linker from PDB
% NH+Nl+NHl = Np; you should check

d_N = 3.8 ;% end-to-end distance in N state along our pulling geometry measured in PDB
%%%%%%%%%%%%%%%%%%%% Inf.of target MP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% contour length per one alpha helix
Lh = zeros(m,1);
Lhl = zeros(m-1,1);
ki = {};
kil = {};
delxh = zeros(m,length(Fmag)); % KR model for each helix
delxl = zeros(m-1,length(Fmag));

for i = 1:m
    Lh(i) = lh*NH(i);
    ki{i} = sqrt(f*(Lh(i)^2)/(4*Ph))';
    for  j = 1:length(Fmag)
        % Kessler-Rabin formula
        delxh(i,j) = (-1./(2*f(j))-ki{i}(j)./(f(j).*tanh(2*ki{i}(j)))+Lh(i)./(tanh(f(j).*Lh(i)))-2*(ki{i}(j).^2).*(1./tanh(f(j)*Lh(i))-f(j).*Lh(i)./(sinh(f(j).*Lh(i))).^2-1)./(3*f(j)));
    end
end


for i = 1:m
    Lhl(i) = lh*NHl(i);
    kil{i} = sqrt(f*(Lhl(i)^2)/(4*Ph))';
    for  j = 1:length(Fmag)
        % Freely Jointed Chain Model -- primary structures in linker part
        %delxl(i,j) = Nl(i).*bk*(coth(f(j).*bk)-1./(f(j).*bk));
        if Lhl(i) == 0
            delxl(i,j) = 0;
        else
            %Kessler-Rabin formula -- secondary structures in linker part
            delxl(i,j) = (-1./(2*f(j))-kil{i}(j)./(f(j).*tanh(2*kil{i}(j)))+Lhl(i)./(tanh(f(j).*Lhl(i)))-2*(kil{i}(j).^2).*(1./tanh(f(j)*Lhl(i))-f(j).*Lhl(i)./(sinh(f(j).*Lhl(i))).^2-1)./(3*f(j)));
        end
    end
end

xl = zeros(length(Fmag),1); % Kessler-Rabin formula (KR model as changing number of helix)
xh = zeros(length(Fmag),1); % Kessler-Rabin formula (KR model as changing number of helix)
hsxGLUT3 = zeros(length(Fmag),1);

xp = zeros(length(Fmag),m+1); % Marko-Siggia fomula (WLC model as changing number of helix)
xhp = zeros(length(Fmag),m); % Kessler-Rabin formula (KR model as changing number of helix)
xGLUT3 = zeros(length(Fmag),m+1);

imGLUT3 = zeros(length(Fmag),1);
% t == number of remainded helices in GLUT3
for j = 1:length(Fmag)
    for tt = 0:m
        % Extensible WLC with streching modulus is approximately infinite
        xp(j,tt+1)= (lp*(Np-sum(NH(1:tt))))*(4/3-4./(3*sqrt(((Fmag(j).*Pp)./kbT)+1))-(10*exp((900./((Fmag(j).*Pp)./kbT)).^(0.25)))./(sqrt(((Fmag(j).*Pp)./kbT))*(exp((900./((Fmag(j).*Pp)./kbT)).^(0.25))-1).^2)+ (((Fmag(j).*Pp)./kbT).^(1.62))./(3.55+3.8*((Fmag(j).*Pp)./kbT).^(2.2))+ Fmag(j)./Kaa );
        % Kessler-Rabin formula for Helix
        xhp(j,tt+1) = sum(delxh(1:tt,j));
        xGLUT3(j,tt+1) = xp(j,tt+1) + xhp(j,tt+1); % consider polypeptide linker + helix mixed state
    end
    
    xl(j) = sum(delxl(1:m-1,j));
    xh(j) = sum(delxh(1:m,j));
    %hsxGLUT3(j) = xl(j)+xh(j);
    hsxGLUT3(j) = xp(j,1)*(sum(Nl))/Np + xl(j) +xh(j); % consider polypeptide linker + helix state out of membrane
    
    imGLUT3(j) = xp(j,1)*(sum(Nl))/Np + xl(j); % consider only linker at 1st stage of membrane protein
    
end

if 0
    figure(1);
    plot(hsxGLUT3-5,Fmag,'m','LineWidth',2);
    hold on;
    plot(xp(:,1)-5,Fmag,'b','LineWidth',2);
    set(gca,'xtick',0:10:125,'ytick',0:5:50,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
    grid on;
    
    xlim([0 125])
    ylim([0 50])
    ylabel('Force (pN)');
    xlabel('Unfolding extension(nm)')
end
