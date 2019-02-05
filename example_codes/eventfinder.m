%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
% Copyright (C): Deyu MING
% Date: 5 Feb. 2019
% Affiliation: Dept of Statistical Science at University College London
% Email: deyu.ming.16@ucl.ac.uk
%
% Reference: Ming, D., Huang, C., Peters, G.W., and Galasso, C. 
% An advanced estimation algorithm for ground-motion models with 
% spatial correlation. BSSA, 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y_event,x_event,w_event]=eventfinder(eventid,id,y,x,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function extracts IMs and covariates of a single event from a list of
%ground-motion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% eventid: is the id that specifies the event to be extracted;
%
% id: a vector of ids that define the events in the ground-motion databse;
%
% y: the IM data at stations of all events in the ground-motion databse; 
%
% x: the covariates at stations of all events in the ground-motion databse;
%
% w: the coordinates (latitude and longtitude) at stations of all events in
% the ground-motion databse.
%
% OUTPUT:
% y_event: the IM data at stations of event labeled by eventid;
%
% x_event: the covariates at stations of event labeled by eventid;
%
% w_event: the coordinates (latitude and longtitude) of stations in
% considered event.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g,event]=findgroups(id);
eventno=find(event==eventid);
indices=find(g==eventno);
y_event=y(indices,:);
x_event=x(indices,:);
w_event=w(indices,:);
end