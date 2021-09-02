%Taking a polynomial estimate of the bottom 10% of the potential surface

% Inputs: 
% PotentialSurface - the potential surface 
% DomainLinspace - the x values for the PotentialSurface
% Degree - the degree of the polynomial fitting, usually 3 
% Proportion - the proportion of the surface that will be fitted 

% Outputs
% coef - the coefficients of the polynomial fitting 

function [coef] = EFM_PolynomialFitting(PotentialSurface,DomainLinspace,Degree,Proportion)

%First determine where the fitting is to be taken 
PS_min = min(PotentialSurface);
PS_max = max(PotentialSurface);
Position_min = find(PotentialSurface==PS_min,1);
Range = PS_max - PS_min;
FittingRange = Proportion*Range; 
FittingHeight = PS_min + FittingRange; %The fitting is taken over the vertical range (PS_min,fitting_height)

%Now determine what values of PotentialSurface are closest to FittingHeight 
PS_before_min = PotentialSurface(1:Position_min);
PS_after_min = PotentialSurface(Position_min:end);
left = find( abs(PS_before_min-FittingHeight)==min(abs(PS_before_min-FittingHeight)),1 );
right = find( abs(PS_after_min-FittingHeight)==min(abs(PS_after_min-FittingHeight)),1 );

%The polynomial fitting will be done on DomainLinspace(LEFT:RIGHT),PotentialSurface(LEFT:RIGHT)
LeftFittingEdge = left;
RightFittingEdge = Position_min + right - 1;

[coef] = polyfit(DomainLinspace(LeftFittingEdge:RightFittingEdge),PotentialSurface(LeftFittingEdge:RightFittingEdge),Degree);

end