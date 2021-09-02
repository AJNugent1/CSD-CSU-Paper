%This code is created for use in the MMath Project: Investigating the
%potetntial of early warning signals in disease elimination. 

%Gillespie simulation of the SIS model 

%{
Inputs: 
I0 - initial number of infected individuals 
N - total population size 
beta - infectivity constant 
gamma - recovery rate 
max_time - time for which the simulation runs 
dt - timestep for incidence aggregation
plots - if 0 then does not plot any graphs, if 1 plots output graphs 

Outputs:
Timeseries for Susceptibles, Infected and Incidence aggregated over a
timestep of length dt 
%}

function [T,I] = GillespieSIS_FixedParameters(para,MaxTime,Startplace)

% Find intial conditions as an integer. 
I0 = Startplace*para.N;

clear I T

%Starting values for S - susceptibles, I - infected 
S(1)=para.N - I0;
I(1)=I0; 

%Initialising time (t) and time vector (T). Index is the position of time t
%in the vector T
T(1)=0;
t=0;
ind=1;

%Main iteration 
while t<MaxTime
    
    %Calculating the total rate of reactions 
    rate(1) = para.Beta*S(ind)*I(ind)/(para.N); %Rate of infections 
    rate(2) = para.Gamma*I(ind); %Rate of recoveries 
    lambda = sum(rate);
    
    if lambda>0 % i.e. there is still happening on in the simulation
        
        %Calculating the timestep 
        r1=rand;
        DT = -log(r1)/lambda;
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
        t=MaxTime;
        T(ind+1)=t;
        S(ind+1)=S(ind);
        I(ind+1)=I(ind);
    end
    
    ind=ind+1;
end

I = I/para.N;

end