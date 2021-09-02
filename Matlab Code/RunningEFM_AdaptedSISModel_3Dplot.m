% Creating the 3D plots of the potential surfaces. This function can only
% be run after the potential surfaces have been calculated for all
% fixed parameter values.

% Inputs:
% n - model control parameter. Expect CSD for n=0,1; CSU for n=2,3. 
% Z_EFM - potential surfaces from the EFM
% Z_Ana - analytic potential surfaces
% Domain - the range over which the potential surfaces are approximated 
% DomainLinspace - the (evenly spaced) points in the domain at which the potential surfaces are
%                  approximated
% NumberSections - the number of (evenly sized) sections the domain is split into
% NumberR0Values - the number of R0 values for which the EFM has been run
% R0_vec - list of R0 values for which the EFM has been run 
% PlotOnFigure - number of the figure on which the 3D plots will be plotted
% OnSubplot - 0: The 3D plots will take up the entire figure
%             1: The 3D plots will be in one subplot on the figure

% Outputs: 
% 3D plots of the potential surfaces 



function RunningEFM_AdaptedSISModel_3Dplot(Z_EFM,Z_Ana,DomainLinspace,Domain,NumberSections,R0_vec,PlotOnFigure,OnSubplot,NumberR0Values,n)

for r = 1:NumberR0Values 

% Controlling axes
ub_3Dplot = max(Z_EFM(end,1:round(NumberSections*0.4)));
lb_3Dplot = min(Z_EFM(1,:));

% Shaping potential surface data
X_3Dplot = [DomainLinspace Domain(2) Domain(1)];
Y1_3Dplot = R0_vec(r)*ones(1,NumberSections);
Y_3Dplot = [Y1_3Dplot R0_vec(r) R0_vec(r)];
Z_3Dplot = [Z_EFM(r,:) lb_3Dplot lb_3Dplot];

% Plotting potential surfaces
figure (PlotOnFigure)
if OnSubplot==1
subplot(2,2,n+1)
end
hold on
f1 = fill3(X_3Dplot,Y_3Dplot,Z_3Dplot,'White');
xb_3Dplot=DomainLinspace( find(Z_EFM(r,:)==min(Z_EFM(r,:)),1) );
yb_3Dplot=R0_vec(r);
zb_3Dplot=min(Z_EFM(r,:));

% Coloring according to R0 value
if R0_vec(r)>1
f1.EdgeColor = 'b';
plot3(xb_3Dplot,yb_3Dplot,zb_3Dplot,'b*')
else 
f1.EdgeColor = [0.8500 0.3250 0.0980];
plot3(xb_3Dplot,yb_3Dplot,zb_3Dplot,'Color',[0.8500 0.3250 0.0980],'Marker','*');
end

% Plotting analytic comparison 
plot3(DomainLinspace,Y1_3Dplot,Z_Ana(r,:),'r--')
xb2_3Dplot=DomainLinspace( find(Z_Ana(r,:)==min(Z_Ana(r,:)),1) );
zb2_3Dplot=min(Z_Ana(r,:));
plot3(xb2_3Dplot,yb_3Dplot,zb2_3Dplot,'r*')

end

% Axes/labels/etc for 3D plot
figure (PlotOnFigure)
if OnSubplot==1
subplot(2,2,n+1)
end
xlabel('Prevalence')
ylabel('R0')
zlabel('Potential Surface')
title(strcat('n=',num2str(n)))
xlim([0 1/2])
zlim([lb_3Dplot ub_3Dplot])
set(gca,'Ydir','reverse')
grid on

end