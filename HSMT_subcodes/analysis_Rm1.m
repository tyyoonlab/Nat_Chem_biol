%% load files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F0 = 0.043; z0 = 25.6; A1 = 35.171; d1 = 0.757; A2 =46.935; d2 =  2.571; % M270 (20191204 for Rm1_after X-actuator setup) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rbead = 1400; T = 296.15;
correction_factor = 0.878; pixel_size = 80;
tmp = dir('*');
pths = arrayfun(@(i) {[tmp(i).name,'/']}, 1:numel(tmp)); npth = numel(pths);
[finfo,fname,nbead,nframe,fps,Roff,ori,f,t,t2,M,R,P,F1,rz,rx,ry,x,y,z,dx,dy,dz] = deal(cell(npth,1));
nfile = zeros(npth,1);

for p = 1 % Bead number
    nfl = 0;
    finfo{p} = dir([pths{p},'r*.xls']);
    nfile(p) = numel(finfo{p});
    [fname{p},Roff{p},ori{p},f{p},M{p},R{p},F1{p},x{p},y{p},z{p},dx{p},dy{p},dz{p}] = deal(cell(nfile(p),1));
    [nbead{p},fps{p},nframe{p}] = deal(zeros(nfile(p),1));
    for n = 1:nfile(p)  % file numver at same bead
        fname{p}{n} = finfo{p}(n).name;
        fname_motor = ['s',fname{p}{n}(2:end)];
        dat = dlmread([pths{p},fname_motor]);
        t2{p}{n} = dat(:,1);
        M{p}{n} = dat(:,2);
        F1{p}{n} = F0 + A1*exp(-(z0-M{p}{n})./d1) + A2*exp(-(z0-M{p}{n})./d2);
        R{p}{n} = (dat(:,3)-floor(dat(:,3)))*360;
        P{p}{n} = dat(:,4);
        
        dat = dlmread([pths{p},fname{p}{n}]);
        nframe{p}(n) = size(dat,1);
        tmp = dlmread([pths{p},'c',fname{p}{n}(2:4),'.fps']);
        fps{p}(n) = tmp(1);
        %         if size(tmp,1)>1
        Roff{p}{n} = tmp(2,:);
        ori{p}{n} = tmp(3,:);
        %         end
        f{p}{n} = dat(:,1);
        dat = dat(:,2:end);
        t{p}{n} = f{p}{n}/fps{p}(n);
        nbead{p}(n) = size(dat,2)/3-1;
        % subtract xy offset
        dat(:,[1:3:end,2:3:end]) = dat(:,[1:3:end,2:3:end]) - repmat(mean(dat(31:60,[1:3:end,2:3:end]),1),[nframe{p}(n),1]);
        rx{p}{n} = dat(:,1)*pixel_size;
        ry{p}{n} = dat(:,2)*pixel_size;
        rz{p}{n} = dat(:,3)*correction_factor;
        x{p}{n} = dat(:,4:3:end)*pixel_size;
        y{p}{n} = dat(:,5:3:end)*pixel_size;
        z{p}{n} = dat(:,6:3:end)*correction_factor;
        dx{p}{n} = (x{p}{n}-repmat(rx{p}{n},[1,nbead{p}(n)]));
        dy{p}{n} = (y{p}{n}-repmat(ry{p}{n},[1,nbead{p}(n)]));
        dz{p}{n} = (z{p}{n}-repmat(rz{p}{n},[1,nbead{p}(n)]));
        
        F1{p}{n} = interp1(t2{p}{n},F1{p}{n},t{p}{n},'linear','extrap');
        R{p}{n} = interp1(t2{p}{n},R{p}{n},t{p}{n},'linear','extrap');
        P{p}{n} = interp1(t2{p}{n},P{p}{n},t{p}{n},'linear','extrap');
        M{p}{n} = interp1(t2{p}{n},M{p}{n},t{p}{n},'linear','extrap');
        % synchronize motor data
    end
    
    clear dat;
    if size(dx{p}{1},2) == 1
        Time{p} = cell2mat(t{p}');
        Tracex{p}(:,1) = cell2mat(dx{p});
        Tracey{p}(:,1) = cell2mat(dy{p});
        Tracez{p}(:,1) = cell2mat(dz{p});
        Force{p} = cell2mat(F1{p});
    else
        sn = size(dx{p}{1},2);
        for sm = 1:sn
            tx = cell2mat(dx{p});
            ty = cell2mat(dy{p});
            tz = cell2mat(dz{p});
            Tracex{p}{sm}(:,1) = tx(:,sm);
            Tracey{p}{sm}(:,1) = ty(:,sm);
            Tracez{p}{sm}(:,1) = tz(:,sm);
        end
        Force{p} = cell2mat(F1{p});
        Time{p} = cell2mat(t{p}');
    end
    Trace{1}(:,1) = Tracex{1};
    Trace{1}(:,2) = Tracey{1};
    Trace{1}(:,3) = Tracez{1};
    disp("File analysis finished. Press any key");
    pause;
end


