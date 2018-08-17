% 
% --> function [A,B,Fe]=Nemoh(w, dir, depth)
%
% Purpose: Matlab wrapper for calculation of hydrodynamic coefficients using Nemoh
% 
% Inputs :
% - w     : Vector length(w) of wave frequencies (rad/s)
% - dir   : Wave direction (degrees)
% - depth : water depth (m), 0 for deep water.
%
% Outputs :
% - A  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of added mass coefficients
% - B  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of radiation damping coefficients
% - Fe : Matrix (6xnBodies)xlength(w) of exciation forces (complex
% values)
%
% Copyright Ecole Centrale de Nantes 2014
% Licensed under the Apache License, Version 2.0
% Written by A. Babarit, LHEEA Lab.
%
function [A,B,Fe]=Nemoh(w,p)
% Preparation du calcul
dir=180;    %head on 
depth=0;    %unlimited depth
name=sprintf('%d',p);   %Individuals in a population are saved in separate folders, repeated each gen
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'ID.dat'));       %change all to \nomrep\nomrep***
line=fgetl(fid);
rep=fscanf(fid,'%s',1);
fclose('all');
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'Nemoh.cal'),'r');%[rep,filesep,'Nemoh.cal']),'r');   %check
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
fclose(fid);
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'Nemoh.cal'),'r'); %check '\' not needed     %[rep,filesep,'Nemoh.cal']),'r');   %check
n=1;
clear textline;
textline={};
while (~feof(fid))
    textline(n)={fgetl(fid)};
    if (n == 4) 
        textline(n)={sprintf('%f                 ! DEPTH			! M		! Water depth',depth)};
    end
    if ((mod(n,18) == 9) && ((n-9)/18 <= nBodies))
        temp=cell2mat(textline(n));
        temp2=[];
        ntemp=length(temp);
        k=1;
        for i=1:ntemp
            if (temp(i) == '\')
                temp2=[temp2,temp(k:i),'\'];
                k=i+1;
            end;            
        end
        temp2=[temp2,temp(k:ntemp)];
        textline(n)={temp2};
        cell2mat(textline(n));
    end
    if (n == 9+18*nBodies) 
        textline(n)={sprintf('%g %f %f       ! Number of wave frequencies, Min, and Max (rad/s)',length(w),w(1),w(length(w)))};
    end
     if (n == 10+18*nBodies) 
        textline(n)={sprintf('%g %f %f		! Number of wave directions, Min and Max (degrees)',1,dir,dir)};
    end
    n=n+1;
end
fclose(fid);
fid = fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'Nemoh.cal'),'w'); %check       %[rep,filesep,'Nemoh.cal']), 'w');    %check
for i=1:n-1
    fprintf(fid, [cell2mat(textline(i)),'\n']);
end
fclose(fid);
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'input.txt'),'wt');%[rep,filesep,'input.txt']),'wt');      %check
fprintf(fid,' \n 0 \n');
status=fclose(fid);
% Calcul des coefficients hydrodynamiques
l = isunix;
if l == 1
    fprintf('\n------ Starting NEMOH ----------- \n');
    system('C:\Users\kmct\Downloads\',name,'preProc');
    fprintf('------ Solving BVPs ------------- \n');
    system('C:\Users\kmct\Downloads\',name,'Solver');
    fprintf('------ Postprocessing results --- \n');
    system('C:\Users\kmct\Downloads\',name,'postProc');
else
    fprintf('\n------ Starting NEMOH ----------- \n');
    text1 = ['C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\Nemoh\preProcessor.exe'];
    system(text1);                          
    fprintf('------ Solving BVPs ------------- \n');
    text2 = ['C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\Nemoh\Solver.exe'];
    system(text2);
    fprintf('------ Postprocessing results --- \n');
    text3 = ['C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\Nemoh\postProcessor.exe'];
    system(text3);
end
%% Lecture des resultats CA CM Fe
clear Periode A B Famp Fphi Fe;
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'Nemoh.cal'),'r');%[rep,filesep,'Nemoh.cal']),'r');   %check
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
for i=1:2+18*nBodies
    ligne=fgetl(fid);
end
nw=fscanf(fid,'%g',1);
fclose(fid);
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'\results','ExcitationForce.tec'),'r');%[rep,filesep,'results',filesep,'ExcitationForce.tec']),'r');   %check
ligne=fgetl(fid);
for c=1:6*nBodies
    ligne=fgetl(fid);
end;
ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+12*nBodies);
    w(i)=ligne(1);
    for j=1:6*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end;
end;
status=fclose(fid);
fid=fopen(fullfile('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',name,'\',sprintf('%d',p),'\results','RadiationCoefficients.tec'),'r');%[rep,filesep,'results',filesep,'RadiationCoefficients.tec']),'r');
ligne=fgetl(fid);
for i=1:6*nBodies
    ligne=fgetl(fid);
end;
for i=1:nBodies*6
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+12*nBodies);
        for j=1:6*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end;
        ligne=fgetl(fid);
    end;
end;
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
i=sqrt(-1);
Fe=Famp.*(cos(Fphi)+i*sin(Fphi));
end