%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%Copyright (C): Deyu MING
%Date: 5 Feb. 2019
%Affiliation: Dept of Statistical Science at University College London
%Email: deyu.ming.16@ucl.ac.uk
%
% Reference: Ming, D., Huang, C., Peters, G.W., and Galasso, C. 
% An advanced estimation algorithm for ground-motion models with 
% spatial correlation. BSSA, 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=shakemap(latRange,lonRang,border,station,predloc,prediction,logscale,lim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function draws the Shake map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% latRange: a vector containing the max and min latitudes;
%
% lonRange: a vector containing the max and min longitudes;
%
% border: a two-column matrix containing the latitude and longitude coordinates
% of the interested region border;
%
% station: a two-column matrix containing the latitude and longitude coordinates
% of stations, the IMs at which are used to estimate GMPE; if station='off',
% the stations are not plotted on the hazard map;
%
% predloc: a two-column matrix containing the latitude and longitude coordinates
% of locations at which the predictions are made;
%
% prediction: the values of IMs at the locations defined in predloc;
%
% logscale: 1: the color map is drawn under logscale; 0: the color map is
% drawn under ordinary scale. Set logscale=1 when you want to plot IM and
% logscale=0 when you want to plot logarithmic IM;
%
% lim: a two-element vector specify the limits of colorbar and colormap.
%
% OUTPUTS:
% Shake map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure;
worldmap(latRange, lonRang);
hold on
if logscale==1
set(gca,'colorscale','log')
end
scatterm(predloc(:,1),predloc(:,2),7,prediction,'s','filled');
colormap('jet')

if ~strcmp(station,'off')
plotm(station(:,1),station(:,2),'^',...
'MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor','k');
end

plotm(border(:,1),border(:,2),'Color',[0.1 0.1 0.1],'linewidth',2);

set(gcf,'units','inches','position',[0,0,5.5,3.5]);
caxis manual
caxis(lim);
cb=colorbar('Location','eastoutside');
cb.Limits=lim;
end
