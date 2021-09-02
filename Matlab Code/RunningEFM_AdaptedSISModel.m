% Main file for running the equation free method on the adapted SIS model. 

% For each combination of n and R0 we approximate the potential surface and
% diffusion function using the equation-free method (EFM). From this we calculate:
% - The location of the fixed point
% - The steepness of the potential surface using a polynomial fitting
% - The value of the diffusion function at the fixed point
% - The value of five common EWS: variance, coefficient of variation,
%   autocorrelation, index of dispersion, decay time.
% All of these values (except the location of the fixed point) are plotted
% against analytic expressions.
% 3D plots of the potential surfaces, again with analytic comparisons, are
% produced for each value of n. 


clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% n-independent setting for the EFM
NumberRealisations = 100;
SimulationMaxTime = 10; % (days)
Timestep = 0.001; % (days)
Domain = [0 0.5];
NumberSections = 50;
RoundStartingValue = 'N'; % Select 'Y' or 'N'
N = 10000; % Population size
DomainLinspace = linspace(Domain(1),Domain(2),NumberSections);
DataSpacing = 'Uniform'; % Select from 'Uniform','Random','RandomDaily'
NumberRandomDataPoints = 100; % Specify if choosing 'Random' for data spacing.
DifferenceMethod = 'Forward'; % Select from 'Forward','Central'

% For each n value we wish to test several R0 values
Number_nValues = 4; % n=0,1,2,3
NumberR0Values = 5;
R0_start = 1.5; % The starting (largest) R0 value
R0_end = 1.1; % The final (smallest) R0 value
R0_vec = linspace(R0_start,R0_end,NumberR0Values); % Vector of R0 values

% Set up empty arrays to store EFM results 
FixedPoint_EFM = zeros(Number_nValues,NumberR0Values);
Steepness_EFM = zeros(Number_nValues,NumberR0Values);
Noise_EFM = zeros(Number_nValues,NumberR0Values);
Variance_EFM = zeros(Number_nValues,NumberR0Values);
CoefVariation_EFM = zeros(Number_nValues,NumberR0Values);
IndexDispersion_EFM = zeros(Number_nValues,NumberR0Values);
Autocorrelation_EFM = zeros(Number_nValues,NumberR0Values);
DecayTime_EFM = zeros(Number_nValues,NumberR0Values);
% Set up empty arrays to store analytic (Ana) results 
FixedPoint_Ana = zeros(Number_nValues,NumberR0Values);
Steepness_Ana = zeros(Number_nValues,NumberR0Values);
Noise_Ana = zeros(Number_nValues,NumberR0Values);
Variance_Ana = zeros(Number_nValues,NumberR0Values);
CoefVariation_Ana = zeros(Number_nValues,NumberR0Values);
IndexDispersion_Ana = zeros(Number_nValues,NumberR0Values);
Autocorrelation_Ana = zeros(Number_nValues,NumberR0Values);
DecayTime_Ana = zeros(Number_nValues,NumberR0Values);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Beginning n loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Begin loop over n values.
% For n=0,1 we expect to see CSD.
% For n=2,3 we expect to see CSU.
for n=0:3
    
    % Since arrays must begin at 0 we use nn for indexing n
    nn = n+1; 
    
    % Create empty arrays to store information for use in 3D plot
    Z_EFM = zeros(NumberR0Values,NumberSections); 
    Z_Ana = zeros(NumberR0Values,NumberSections); 
    Beta_vec = zeros(1,NumberR0Values);
    Gamma_vec = zeros(1,NumberR0Values);
    
    % Begin loop over R0 values. 
    % THIS IS A PARFOR LOOP. This can be changed to a normal 'for' loop. 
    
    parfor r=1:NumberR0Values
        
        % Set R0
        R0 = R0_vec(r);
        
        % Set Beta and Gamma:
        
        % First calculate the normalisation constand K_n
        K = ( R0_start/(R0_start - 1) )^(-n);
        % Then calculate beta and gamma
        Gamma = K*( (R0/(R0-1))^n );
        Beta = R0*Gamma;
        % Store these values to use in the 3D plot
        Beta_vec(r) = Beta;
        Gamma_vec(r) = Gamma;
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running EFM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Set up the EFM
        SimulationParameters = struct('N',N,'Beta',Beta,'Gamma',Gamma);
        SimulationFunction = @(SimulationParameters,SimulationMaxTime,StartingValue)...
                            GillespieSIS_FixedParameters(SimulationParameters,SimulationMaxTime,StartingValue);
        
                        
        % Run the EFM
        [DriftFunction_EFM,DiffusionFunction_EFM] = EquationFreeMethod(SimulationFunction,...
                                                                       SimulationParameters,...
                                                                       SimulationMaxTime,...
                                                                       NumberRealisations,...
                                                                       Timestep,...
                                                                       Domain,...
                                                                       RoundStartingValue,...
                                                                       NumberSections,...
                                                                       DataSpacing,...
                                                                       NumberRandomDataPoints,...
                                                                       DifferenceMethod);
        
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations from EFM %%%%%%%%%%%%%%%%%%%%%%%%%%
                                       
                                       
        % We now calculate the potential surface, fixed point location,
        % steepness of potential surface and diffusion at the fixed point.
                                       
        % Integrate -F_est numerically to give potential surface
        PotentialSurface_EFM = cumtrapz(DomainLinspace,-DriftFunction_EFM);
        Z_EFM(r,:) = PotentialSurface_EFM; % For 3D plot
        
        % Find the location of the fixed point
        % First find the minimum of the potential surface
        PositionMinimum = find(Z_EFM(r,:)==min(Z_EFM(r,:)),1);
        % Then the corresponding value in the domain
        FixedPoint_EFM(nn,r) = DomainLinspace( PositionMinimum );
        
        % Calculate polynomial (cubic) fitting to bottom of potential surface
        [coef] = EFM_PolynomialFitting(PotentialSurface_EFM,DomainLinspace,3,0.1);
        % Then calculate the steepness of the cubic at the fixed point
        Steepness_EFM(nn,r) = 6*coef(1)*FixedPoint_EFM(nn,r) + 2*coef(2);

        % Calculate the value of the diffusion function at fixed point. 
        % (This array is called Noise_EFM to avoid confusion with the
        % diffusion function itself) 
        Noise_EFM(nn,r) = DiffusionFunction_EFM(PositionMinimum);
        
        % We also calcualte the analytic potential surface.
        DriftFunction_Ana = Beta.*DomainLinspace.*(1 - DomainLinspace) - Gamma.*DomainLinspace;
        PotentialSurface_Ana = cumtrapz(DomainLinspace,-DriftFunction_Ana);
        Z_Ana(r,:) = PotentialSurface_Ana; % For 3D plot
        
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Creating 3D plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot the 3D potential surfaces figure. This is done within the
    % n-values for loop. 
    
    % Specify which figure the 3D plots should appear on. 
    PlotOnFigure = 1; 
    OnSubplot = 1;
    % Separate function to plot 3D plots
    RunningEFM_AdaptedSISModel_3Dplot(Z_EFM,Z_Ana,DomainLinspace,Domain,...
                                      NumberSections,R0_vec,PlotOnFigure,OnSubplot,NumberR0Values,n)
    
                                  
%%%%%%%%%%%%%%%%%%%%%%%% Calculating EWS from EFM %%%%%%%%%%%%%%%%%%%%%%%%%


    Variance_EFM(nn,:) = (Noise_EFM(nn,:).^2)./(2.*Steepness_EFM(nn,:));
    CoefVariation_EFM(nn,:) = Noise_EFM(nn,:)./(FixedPoint_EFM(nn,:).*sqrt(2.*Steepness_EFM(nn,:)));
    IndexDispersion_EFM(nn,:) = (Noise_EFM(nn,:).^2)./(2.*Steepness_EFM(nn,:).*FixedPoint_EFM(nn,:));
    Autocorrelation_EFM(nn,:) = exp(-Steepness_EFM(nn,:));
    DecayTime_EFM(nn,:) = 1./Steepness_EFM(nn,:);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%% Analytic calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%
                      

    FixedPoint_Ana(nn,:) = 1 - 1./R0_vec;
    K = ( R0_start/(R0_start - 1) )^(-n);
    Steepness_Ana(nn,:) = (R0_vec.*K).*((R0_vec./(R0_vec - 1)).^(n-1));
    Noise_Ana(nn,:) = sqrt((2./(N/K)).*((R0_vec./(R0_vec - 1)).^(n-1)));
    Variance_Ana(nn,:) = (Noise_Ana(nn,:).^2)./(2.*Steepness_Ana(nn,:));
    CoefVariation_Ana(nn,:) = Noise_Ana(nn,:)./(FixedPoint_Ana(nn,:).*sqrt(2.*Steepness_Ana(nn,:)));
    IndexDispersion_Ana(nn,:) = (Noise_Ana(nn,:).^2)./(2.*Steepness_Ana(nn,:).*FixedPoint_Ana(nn,:));
    Autocorrelation_Ana(nn,:) = exp(-Steepness_Ana(nn,:));
    DecayTime_Ana(nn,:) = 1./Steepness_Ana(nn,:);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%% Finished this n value %%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % Display that this n value is complete.
    disp(strcat('Completed n=',num2str(n)))
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Provides consistent colouring for n=0,1,2,3
colours = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];

LineStyleEFM = '-';
MarkerStyleEFM = 'none';
LineWidthEFM = 0.5;
LineStyleAna = '--';
MarkerStyleAna = 'none';
LineWidthAna = 0.5;

% Steepness
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0_vec,Steepness_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    plot(R0_vec,Steepness_Ana(nn,:),'Color',colours(nn,:),'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Steepness of potential surface')
legend('Location','northwest','FontSize',12)

% Noise
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0_vec,Noise_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    plot(R0_vec,Noise_Ana(nn,:),'Color',colours(nn,:),'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Diffusion at fixed point')
legend('Location','northwest','FontSize',12)

% Variance
figure
hold on
for n=0:3
    nn=n+1;
    plot(R0_vec,Variance_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    if n==3
    plot(R0_vec,Variance_Ana(nn,:),'Color',[0 0 0],'DisplayName',strcat('Analytic, n=',num2str(n)),...
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
    plot(R0_vec,CoefVariation_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    if n==3
    plot(R0_vec,CoefVariation_Ana(nn,:),'Color',[0 0 0],'DisplayName',strcat('Analytic, n=',num2str(n)),...
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
    plot(R0_vec,IndexDispersion_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    if n==3
    plot(R0_vec,IndexDispersion_Ana(nn,:),'Color',[0 0 0],'DisplayName',strcat('Analytic, n=',num2str(n)),...
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
    plot(R0_vec,Autocorrelation_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    plot(R0_vec,Autocorrelation_Ana(nn,:),'Color',colours(nn,:),'DisplayName',strcat('Analytic, n=',num2str(n)),...
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
    plot(R0_vec,DecayTime_EFM(nn,:),'Color',colours(nn,:),'DisplayName',strcat('EFM, n=',num2str(n)),...
        'LineStyle',LineStyleEFM,'Marker',MarkerStyleEFM,'LineWidth',LineWidthEFM)
    plot(R0_vec,DecayTime_Ana(nn,:),'Color',colours(nn,:),'DisplayName',strcat('Analytic, n=',num2str(n)),...
        'LineStyle',LineStyleAna,'Marker',MarkerStyleAna,'LineWidth',LineWidthAna)
end
set(gca,'Xdir','reverse')
xlabel('R0')
ylabel('Decay time')
legend('Location','northwest','FontSize',12)