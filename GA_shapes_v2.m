% GA_shapes_v2.m     % Genetic algorithm for optimal miniWEC float shape KM

clear all; close all; clc

%% GA parameters 
Nind = 30;                          % micro GA w/ 30 individuals
NVAR = 26;                          % number of variables
PRECI = 5;                          % discrete increments - 32 with 0.0098m spacing 
MAXGEN = 25;%50;                        % increase MAXGEN after preliminary runs
GGAP = 0.9;                         % generational gap 

% field descriptor
FieldD = [rep([PRECI],[1,NVAR]);rep([0.9144;1.2192],[1,NVAR]);...
    rep([1;0;1;1],[1,NVAR])];                               % Gray coded

% initial population
Chrom = crtbp(Nind,NVAR*PRECI);

% seed initial population with ellipse, rev disk, and rev conical
Chrom(1,:) = [0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,1,0,1,0,1,1,0,0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,1,0,1,1,0,1,0,0,1,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0,0,1,0,1,0,1,1,0,1,0,1,1,1,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0]; %ellipse
Chrom(2,:) = [1,0,0,0,0,1,0,0,1,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,0,1,1,0,1,1,0,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,1,1,0,0,0,1,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,1,1,1,1,1,0,1,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,0,0,0,0]; %reverse disk
Chrom(3,:) = [1,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0,1,1,0,1,0,1,1,1,1,0,1,0,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,0,0,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0]; %reverse conical 

% Reset counters
   Best = NaN*ones(MAXGEN,1);	% best in current population
   global R_data P_data Z_data
   gen = 0;                     % generational counter
   R_data = zeros(1,28);        % stores radius data 
   Z_data = zeros(1,28);
   P_data = 0;                  % stores power proxy data 
   go = char('go');             % read in 'stop' or 'go' within loop
   tic;

parpool parallel computing
   parpool('local',4);                                         %parallel computing (cloud)
   options = optimoptions('ga','UseParallel', true, 'UseVectorized', false);    %parallel computing GA
    
% Evaluate initial population, call NEMOH in objfunShapes.m 
   ObjV = objfunShapes2(bs2rv(Chrom,FieldD));  

% Generational Loop 
while (gen < MAXGEN & char(go) == 'go')
    
    % read text file -> 'stop' or 'go' to terminate in real time
    [go] = textread('StopGA.dat',...
        '%s');
    
    % Assign fitness values to population (minimizing based on '50000-Proxy.*10')
    FitnV = ranking(ObjV);          
    
    % Select individuals for breeding (stochastic universal sampling)
    SelCh = select('sus', Chrom, FitnV, GGAP);       
    
    % Recombine individuals (single point crossover) 
    SelCh = recombin('xovsp',SelCh,0.7);     
    
    % Mutate
    SelCh = mut(SelCh);  
    
    % Evaluate Offspring, call NEMOH in objfunShapes.m
    ObjVSel = objfunShapes2(bs2rv(SelCh,FieldD));
    
    % Reinsert Offspring into population 
    [Chrom, ObjV] = reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel);   
    
    % Plot best shape fitness function & profile
    clear prof
    [Best(gen+1),I] = min(ObjV);
    figure(1)
    subplot(1,2,1)
    plot(Best,'bo'); 
    hold on 
    plot(49423,'rx')        %best Proxy of preliminary shapes - teardrop %check**
    xlabel('generation'); ylabel('fitfun(x)');
    text(0.5,0.95,['Best = ', num2str(round(Best(gen+1)))],'Units','normalized'); 
    hold off
    subplot(1,2,2)
    prof = bs2rv(Chrom(I,:),FieldD);   %plot best shape of each generation   
    area(prof);    
   
    % Increment counter
    gen = gen+1;
    
end

toc 