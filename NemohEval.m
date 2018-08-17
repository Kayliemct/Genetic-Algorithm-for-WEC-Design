% function NemohEval.m - NEMOH with parallel computing %KM

function Proxy = NemohEval(p)

    datadir = ['C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',sprintf('%d',p)];
    cd (datadir);        %change directory for each individual
        
% evaluate using NEMOH - NemohGA.m
    w = linspace((2*pi()./8),(2*pi()./0.5));    
    clear A B Fe
    [A,B,Fe] = NemohGA2(w,p);           

% get desired added mass and damping data 
    for jj=1:6
        a(jj,:) = A(jj,jj,:);
        b(jj,:) = B(jj,jj,:);
    end

% Filter NEMOH data - heave only 
    af = medfilt1(a(3,:));     
    bf = medfilt1(b(3,:));
    FefA = medfilt1(abs(Fe(:,3)));
    Fefph = medfilt1(angle(Fe(:,3)))+pi();

% Power proxy = a(MaxSj)+b(DxSj)+c(ExxSj)  
    Hmo = 0.18;        
    Tp = 1.5;           
    gamma = 3.3;        %3.3 for Jonswap
    alpha = 5.061*(Hmo^2/Tp^4)*(1-0.287*log(gamma));
    Sj = ((alpha.*9.81.^2)./(w.^5)).*exp(-(5/4).*(3.5./w).^4).*gamma.^exp((-(w./3.5-1).^2)./(2.*0.07.^2));
    T = (2.*pi())./w;

    M = trapz((af.*Sj)./T.^2);
    D = trapz(bf.*Sj.*T);
    E = trapz(FefA.*Sj'.*sin(Fefph));
   
    coeff = [0.4; -0.1; 0.2];
    Proxy = M.*coeff(1)+D.*coeff(2)+E.*coeff(3);  

end