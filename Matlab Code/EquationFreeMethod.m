% Code for performing the equation-free method on timeseries data

% Inputs: 
% SimulationFunction - the simulation to be performed
% SimulationParameters - the parameters for the simulation
% NumberRealisations - the number of simulations to be run
% Timestep - desired timestep between data points 
% Domain - the domain over which the drift & diffusion functions are
%          approximated. Should be a two column array giving the start and
%          end points of an interval domain. 
% IntegerStartingValue - 0: the starting value does not need to be rounded
%                      - 1: the starting value should be rounded to the
%                           nearest integer
% NumberSections - the number of sections the Domain is split into


% Requirements on SimulationFunction:
% Must take as inputs simulation parameters and a starting position.
% Must output a timeseries in the form (T,X) where T is the timepoints and
% X gives the value of the variable at each timepoint.


% Outputs: 
% F_est - an approximation to the drift function. Gives the approximate
%         value of the drift function in each Section.
% D_est - an approximation to the diffusion function. Gives the approximate
%         value of the diffusion function in each Section.


function [F_est,D_est] = EquationFreeMethod(SimulationFunction,...
                                            SimulationParameters,...
                                            SimulationMaxTime,...
                                            NumberRealisations,...
                                            Timestep,...
                                            Domain,...
                                            RoundStartingValue,...
                                            NumberSections,...
                                            DataSpacing,...
                                            NumberRandomDataPoints,...
                                            DifferenceMethod)

                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            

% Setting up empty arrays
topf=zeros(1,NumberSections);
topd=zeros(1,NumberSections);
bot=zeros(1,NumberSections);
F_est=zeros(1,NumberSections);
D_est_step1=zeros(1,NumberSections);
D_est=zeros(1,NumberSections);


% Setting some notation
M = Domain(2) - Domain(1);
SectionWidth = M/NumberSections;


%%%%%%%%%%%%%%%%%%%%%% Main loop over realisations %%%%%%%%%%%%%%%%%%%%%%%%


for rep=1:NumberRealisations
    
    
%%%%%%%%%%%%%%%%%%%%%%%% Running each realisation %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Set a random starting point within the domain
    StartingValue=M*rand; 
    % Round to nearest integer if needed (e.g. for Gillespie simulations)
    if strcmp('Y',RoundStartingValue)==1
        StartingValue = round(StartingValue);
    end
    
    % Run the simulation
    [T,X]=SimulationFunction(SimulationParameters,SimulationMaxTime,StartingValue); 
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spacing the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if strcmp('Uniform',DataSpacing)==1
        % Linearly space data points 
        TSpaced = linspace(0,SimulationMaxTime,SimulationMaxTime/Timestep);
        XSpaced = interp1(T,X,TSpaced,'previous');
    end
    
    if strcmp('Random',DataSpacing)==1
        if length(T)>=NumberRandomDataPoints
            Indexes = sort(randperm(length(T),NumberRandomDataPoints),'ascend');
        else
            Indexes = sort(randperm(length(T),length(T)),'ascend');
        end
        TSpaced = T(Indexes);
        XSpaced = X(Indexes);
        Timestep = mean(diff(TSpaced));
    end
    
    if strcmp('RandomDaily',DataSpacing)==1
        TDaily = linspace(0,SimulationMaxTime-1,SimulationMaxTime);
        RandomNumbers = rand(1,SimulationMaxTime);
        TSpaced = TDaily + RandomNumbers;
        XSpaced = interp1(T,X,TSpaced,'previous');
        Timestep =1;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Adding to EFM sums %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % We work through X one timepoint at a time. The method depends on the
    % choice of DifferenceMethod.
    
    if strcmp('Forward',DifferenceMethod)==1
        for t=1:length(XSpaced) - 1
            % Determine which Section X lies in at this timepoint 
            for s=1:NumberSections
                if XSpaced(t)>=Domain(1) + (s-1)*SectionWidth && XSpaced(t)<=Domain(1) + s*SectionWidth 
                    % Add to the appropriate sums
                    topf(s)=topf(s) + ( XSpaced(t+1) - XSpaced(t) )/(TSpaced(t+1)-TSpaced(t)); 
                    topd(s)=topd(s) + ( (XSpaced(t+1) - XSpaced(t))^2 )/(TSpaced(t+1)-TSpaced(t));
                    % Add one to the denomiator of the sum
                    bot(s)=bot(s)+1;
                end
            end
        end
    end
    
    if strcmp('Central',DifferenceMethod)==1
        for t=2:length(XSpaced) - 1
            % Determine which Section X lies in at this timepoint 
            for s=1:NumberSections
                if XSpaced(t)>=Domain(1) + (s-1)*SectionWidth && XSpaced(t)<=Domain(1) + s*SectionWidth 
                    % Add to the appropriate sums
                    topf(s)=topf(s) + ( XSpaced(t+1) - XSpaced(t-1) )/(TSpaced(t+1)-TSpaced(t-1)); 
                    topd(s)=topd(s) + ( (XSpaced(t+1) - XSpaced(t-1))^2 )/(TSpaced(t+1)-TSpaced(t-1));
                    % Add one to the denomiator of the sum
                    bot(s)=bot(s)+1;
                end
            end
        end
    end

end %End of simulations 


%%%%%%%%%%%%%%%%%%%%%%%%%%% Averaging EFM sums %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for s=1:NumberSections
    
    if bot(s)>0 
    F_est(s)=topf(s)/bot(s);
    D_est_step1(s) = topd(s)/bot(s);
    else
        if s>1
        F_est(s)=F_est(s-1);
        D_est_step1(s) = D_est_step1(s-1);
        else 
            F_est(s)=0;
            D_est_step1(s)=0;
        end
    end  
    
end


% Calculate the approximation to the diffusion function using F_est as an
% adjustment. 
if strcmp('Forward',DifferenceMethod)==1
    D_est = sqrt( D_est_step1 - (F_est.^2).*Timestep );
end
if strcmp('Central',DifferenceMethod)==1
    D_est = sqrt( D_est_step1 - (F_est.^2).*Timestep.*2 );
end

end