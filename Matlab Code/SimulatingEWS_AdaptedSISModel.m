% Verifying the behaviour of 5 early warning signals for the adapted SIS 
% model through time-varying parameter simulations. 

% EWS are calculated between realisations of Gillespie simulations in which 
% the value of R0 is slowly reduced. These results are compared to analytic
% expressions for the EWS. 


clear 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% n-independent settings and parameters
N = 10000;                      % Population size
R0_start = 1.5;                 % Starting value of R0
R0_end = 1.1;                   % Final value of R0
I0=round(N*(1 - 1/R0_start));   % initial level of infection (at steady state).

NumberRealisations = 100; % The number of realisations of the simulation for each value of n
MaxTime=200;              % Length of each simulation
BurnTime=5;               % Burn-in time before R0 begins to change
RelativeNumberPoints = 1; % The number of data points taken per day.
                          % Equivalently a data point is taken every
                          % 1/RelativeNumberPoints days.
NumberTimePoints = MaxTime*RelativeNumberPoints; 
Number_nValues = 4;       % We test n=0,1,2,3


% Set up empty arrays to store simulated (Sim) results 
Variance_Sim = zeros(Number_nValues,NumberTimePoints);
CoefVariation_Sim = zeros(Number_nValues,NumberTimePoints);
IndexDispersion_Sim = zeros(Number_nValues,NumberTimePoints);
Autocorrelation_Sim = zeros(Number_nValues,NumberTimePoints);
DecayTime_Sim = zeros(Number_nValues,NumberTimePoints);
FixedPoint_Sim = zeros(Number_nValues,NumberTimePoints);
% Set up empty arrays to store analytic (Ana) results 
Variance_Ana = zeros(Number_nValues,NumberTimePoints);
CoefVariation_Ana = zeros(Number_nValues,NumberTimePoints);
IndexDispersion_Ana = zeros(Number_nValues,NumberTimePoints);
Autocorrelation_Ana = zeros(Number_nValues,NumberTimePoints);
DecayTime_Ana = zeros(Number_nValues,NumberTimePoints);
FixedPoint_Ana = zeros(Number_nValues,NumberTimePoints);

% Create a vector of R0 values that matches how R0 decreases during
% simulations (including the burn in time).
NumberBurnInPoints = round(RelativeNumberPoints*BurnTime)-1;
R0Lin(1:NumberBurnInPoints) = R0_start;
R0Lin(NumberBurnInPoints+1:NumberTimePoints) = linspace(R0_start,R0_end,NumberTimePoints - NumberBurnInPoints);


%%%%%%%%%%%%%%%%%%%%%%%%%% Main simulations loop %%%%%%%%%%%%%%%%%%%%%%%%%%


% Begin loop over n values.
% For n=0,1 we expect to see CSD.
% For n=2,3 we expect to see CSU.
for n=0:3
    
    % Since arrays must begin at 1 but n begins at 0 we use nn for indexing
    nn = n+1;
    
    % n-dependent parameters
    k = (R0_start/(R0_start-1))^n;     % normalisation constant
    Gamma = @(r) (1/k).*(r./(r-1)).^n; % gamma rate function
    Beta = @(r) r.*Gamma(r);           % beta rate function
    Parameters = struct('N',N,'Beta',Beta,'Gamma',Gamma,'I0',I0,'MaxTime',MaxTime,...
                        'BurnTime',BurnTime,'R0_start',R0_start,'R0_end',R0_end);
    
    % Creating empty arrays to store data
    IStore = zeros(NumberTimePoints,NumberRealisations);   %empty array to store infection data
    TLin =linspace(0,MaxTime,NumberTimePoints);            %for linearly spacing data
    ILin = zeros(size(TLin));                              %for storing linearly spaced data
    
    % Begin loop over realisations
    % THIS IS A PARFOR LOOP (can be changed to a 'for' loop, but will be slower) 
    parfor Realisation=1:NumberRealisations
        
        % Run the simulation. Produces results for disease prevalence (I)
        % as a PROPORTION of the population. 
        [T,I] = GillespieSIS_TimeVaryingParameters(Parameters);
        
        % Linearly space data points 
        ILin = interp1(T,I,TLin,'previous');
        
        %Store data
        IStore(:,Realisation) = ILin;
        
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%% Calculate EWS - Analytic %%%%%%%%%%%%%%%%%%%%%%%%%


    % Create a vector of Beta and Gamma values that matches how they change 
    % during simulations (including the burn in time).
    BetaLin = Beta(R0Lin);
    GammaLin = Gamma(R0Lin);
    
    % Variance
    Variance_Ana(nn,:) = 1./(N*R0Lin);

    % Lag-1 autocorrelation
    Autocorrelation_Ana(nn,:) = exp(-abs(BetaLin - GammaLin));

    % Decay time
    DecayTime_Ana(nn,:) = (k./R0Lin).*( (R0Lin./(R0Lin-1)).^(1-n) );

    % Location of fixed point
    FixedPoint_Ana(nn,:) = 1 - 1./R0Lin;

    % Coefficient of variation
    CoefVariation_Ana(nn,:) = sqrt(R0Lin./N)./(R0Lin - 1);

    % Index of dispersion
    IndexDispersion_Ana(nn,:) = 1./((R0Lin - 1).*N);
    
    
%%%%%%%%%%%%%%%%%%%%%%%% Calculate EWS - Simulated %%%%%%%%%%%%%%%%%%%%%%%%


    % The mean gives the location of the fixed point
    FixedPoint_Sim(nn,:) = mean(IStore,2);
    
    % Detrend data by taking away the mean between simulations 
    Detrended_Sim = IStore - mean(IStore,2);
    
    % Variance
    Variance_Sim(nn,:) = var(Detrended_Sim,0,2);
    
    % Lag-1 autocorrelation
    for j=3:length(TLin)
        CorrcoefMatrix = corrcoef(IStore(j,:),IStore(j-1,:));  
        Autocorrelation_Sim(nn,j) = CorrcoefMatrix(1,2);
    end
    Autocorrelation_Sim(nn,1) = Autocorrelation_Sim(nn,3);
    Autocorrelation_Sim(nn,2) = Autocorrelation_Sim(nn,3);
    
    % Decay Time
    DecayTime_Sim(nn,:) = -1./log(min(max(Autocorrelation_Sim(nn,:), 0),1));
    
    % Coefficient of Variation
    Standard_Deviation = sqrt(Variance_Sim(nn,:));
    CoefVariation_Sim(nn,:) = Standard_Deviation./FixedPoint_Sim(nn,:);
    
    % Index of Dispersion
    IndexDispersion_Sim(nn,:) = Variance_Sim(nn,:)./FixedPoint_Sim(nn,:);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%% Finished this n value %%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % Display that this n value is complete.
    disp(strcat('Completed n=',num2str(n)))
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Provides consistent colouring for n=0,1,2,3
colours = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];

LineStyleSim = '-';
MarkerStyleSim = 'none';
LineWidthSim = 0.5;
LineStyleAna = '--';
MarkerStyleAna = 'none';
LineWidthAna = 0.5;

% Variance
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0Lin,Variance_Sim(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleSim,'Marker',MarkerStyleSim,'LineWidth',LineWidthSim)
    if n==3
    plot(R0Lin,Variance_Ana(nn,:),'Color',[0 0 0],'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
    end
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Variance')
legend('Location','northwest','FontSize',12)

% Coefficient of variation
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0Lin,CoefVariation_Sim(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleSim,'Marker',MarkerStyleSim,'LineWidth',LineWidthSim)
    if n==3
    plot(R0Lin,CoefVariation_Ana(nn,:),'Color',[0 0 0],'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
    end
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Coefficient of variation')
legend('Location','northwest','FontSize',12)

% Index of dispersion
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0Lin,IndexDispersion_Sim(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleSim,'Marker',MarkerStyleSim,'LineWidth',LineWidthSim)
    if n==3
    plot(R0Lin,IndexDispersion_Ana(nn,:),'Color',[0 0 0],'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
    end
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Index of dispersion')
legend('Location','northwest','FontSize',12)

% Autocorrelation
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0Lin,Autocorrelation_Sim(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleSim,'Marker',MarkerStyleSim,'LineWidth',LineWidthSim)
    plot(R0Lin,Autocorrelation_Ana(nn,:),'Color',colours(nn,:),'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Lag-1 autocorrelation')
legend('Location','southwest','FontSize',10)

% Decay time
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0Lin,DecayTime_Sim(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleSim,'Marker',MarkerStyleSim,'LineWidth',LineWidthSim)
    plot(R0Lin,DecayTime_Ana(nn,:),'Color',colours(nn,:),'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Decay time')
legend('Location','northwest','FontSize',12)