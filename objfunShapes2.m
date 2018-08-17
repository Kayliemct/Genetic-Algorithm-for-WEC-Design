% objfunShapes.m - objective function to use with GA_shapes_v1.m KM

% check if shape has already been evaluated & locate proxy 
% run NEMOH using 'axiMeshGA.m' & 'NemohGA.m' 
% low pass filter on data 
% power proxy = a(MaSj/T^2)+b(DSjT)+c(a(Ex)Sj*sin(phase)) 
% fitness function = abs(50000-proxy.*10)'

function ObjValS = objfunShapes2(Chrom,FieldD);

global R_data P_data Z_data

% Compute population parameters
   [Nind,Nvar] = size(Chrom);
%    save('Chrom.mat','Chrom');    %overwrites previous population chrom.mat
   
% pre-NEMOH evaluation
for i = 1:Nind
   
    datadir = ['C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\',sprintf('%d',i)];
    cd (datadir);        %change directory for each individual

    % radius data points
    r = [0, Chrom(i,1:26), 0];  
        
% vertical spacing for float & number of points for NEMOH    
    z = [0.8128, linspace(0.8128,0,26), 0];
    n = 28;

% find draft for equivalent submerged volume through numerical integration
    Vol_Z = 0;              % for summing volume of the slices
    dZ = z(1,end-4)-z(1,end-3);     % height of slices
    rad = fliplr(r);      % start num. int. from base of shape

    for k=2:length(z)-2
    Radius_k = rad(1,k);   
    Area_k = pi().*Radius_k.^2;
    Volume_k = Area_k.*dZ;
    Vol_Z = Vol_Z + Volume_k;
    if Vol_Z >2.145                % select the draft closest to volume 2.145 m^3
        error1 = Vol_Z-2.145;
        error2 = Vol_Z-Volume_k-2.145;
        if abs(error1) < abs(error2)
            d = k.*dZ;
            break
        else
            d = (k-1).*dZ;
            break
        end
    end
    end

% subtract draft from vertical spacing to get z & find vertical centroid
    z = z-d;
    polyin=polyshape({r},{z});
    [x_c,zCG]=centroid(polyin);        
    
% evaluate using NEMOH pre-processing - axiMeshGA.m
    [m,M_nemoh,KH,XB,YB,ZB] = axiMeshGA2(r,z,n,zCG,i);   
    
% saving z data for each shape
    u0 = size(R_data,1);
    u = u0+i+1;     %space between generations
    Z_data(u,:) = z;
end   

% back to main directory
    cd ('C:\Users\labuser\Documents\Kaylie\GA\100GenParallel\');

% parallel loop running NEMOH 
parfor p = 1:Nind  
    
    Proxy(p) = NemohEval(p); 
    
end

% Fitness function to be minimized  
    ObjValS = abs(50000-Proxy.*10)'; %max Proxy = 5,000 with this fitness function
    
% save r and Proxy for future repeated shapes to solve faster 
    u0 = size(R_data,1);
    for j=1:Nind
        u = u0+j+1;                 % space between offspring generations
        R_data(u,:) = [0, Chrom(j,:), 0];
        P_data(u) = Proxy(j);
    end
end