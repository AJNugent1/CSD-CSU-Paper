%Gillespie simulation of the SIS model with time-varying parameters

% Inputs: 
% I0 - initial number of infected individuals 
% N - total population size 
% beta - infectivity function, this is a function of R0 
% gamma - recovery function, this is a function of R0 
% max_time - time for which the simulation runs 
% burn_time - time for which R0 is constant 
% R00 - initial value of R0
% R01 - final value of R0
% dt - timestep for incidence aggregation
% plots - if 0 then does not plot any graphs, if 1 plots output graphs 

% Outputs:
% Timeseries for the proportin of infected individuals. 

function [T,I] = GillespieSIS_TimeVaryingParameters(para)

clear I T

%Starting values for S - susceptibles, I - infected 
S(1)=para.N - para.I0;
I(1)=para.I0; 

%Initialising time (t) and time vector (T). Index is the position of time t
%in the vector T
T(1)=0;
t=0;
ind=1;

%Main iteration 
while t<para.MaxTime
    
    %Calculating the total rate of reactions 
    
    if t<para.BurnTime
        R0 = para.R0_start;
    else
        R0 = para.R0_start - (para.R0_start-para.R0_end)*(t-para.BurnTime)/(para.MaxTime-para.BurnTime);
    end
    
    b = para.Beta(R0);
    g = para.Gamma(R0);
    
    rate(1) = b*S(ind)*I(ind)/(para.N); %Rate of infections 
    rate(2) = g*I(ind); %Rate of recoveries 
    lambda = sum(rate);
    
    if lambda>0 % i.e. there is still happening on in the simulation
        
        %Calculating the timestep 
        r1=rand;
        DT = -log(r1)/lambda;
        DT = max(DT,1e-08);
        t = t+DT;
        T(ind+1)=t;
        
        %Calculating the next event 
        r2 = rand*lambda;
        event = find(cumsum(rate)>=r2,1);
        if event == 1 %the next event is an infection
            S(ind+1)=S(ind)-1;
            I(ind+1)=I(ind)+1;   
        else 
            if event == 2 % the next event is a recovery
            S(ind+1)=S(ind)+1;
            I(ind+1)=I(ind)-1;
            end
        end
        
    else % i.e. there is nothing going on in the simulation, the disease has died out
        t=para.MaxTime;
        T(ind+1)=t;
        S(ind+1)=S(ind);
        I(ind+1)=I(ind);
    end
    
    ind=ind+1;
end

% Divide by N to give a timeseries for the proportion of infected
% individuals. 
I = I/para.N;

end