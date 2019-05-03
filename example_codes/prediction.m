%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%Copyright (C): Deyu MING
%Date: 3 May 2019
%Affiliation: Dept of Statistical Science at University College London
%Email: deyu.ming.16@ucl.ac.uk
%
% Reference: Ming, D., Huang, C., Peters, G.W., and Galasso, C. 
% An advanced estimation algorithm for ground-motion models with 
% spatial correlation. BSSA, 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z_hat=prediction(y_event,x_event,u_event,w_event,v_event,output,g,cf,base)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the linear predictions at grid points given the
%observations at stations during a specific event and the estimated
%parameters for the model proposed by Akkar and Bommer (2010).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% y_event: the IM data at stations of the event under the consideration;
%
% x_event: the covariates at stations of event under the consideration;
%
% u_event: the covariates at grid points in event under the consideration;
%
% w_event: the coordinates (latitude and longtitude) of stations in
% considered event;
%
% v_event: the coordinates (latitude and longtitude) of grid points in
% considered event;
%
% output: the output structure produced by the Scoring estimation approach;
%
% g: a function handle of the user-defined design matrix function for the
% linear coefficients;
%
% cf: the type of the covariance function: currently support no spatial correlation
% ('No'), exponential type ('Exp'), Matern type with v=1.5 ('Matern1.5') 
% and squared exponential type ('SExp');
%
% base: the base ('10' or 'exp') used to compute the logarithmic IMs. For example, in Akkar
% and Bommer (2010), the log(PGA) is to the base of 10. Set base to 'off' so
% the predictions are produced with log in regardless of the base.
%
% OUTPUT:
% z_hat: predictions at grid points in the event under the consideration
% given the observations at stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=output.ParameterEstimates.beta;
gamma=output.ParaEstimates.gamma;
theta=output.ParaEstimates.theta;
%compute the means at stations in the considered event given estimates of b
fx_event=g(x_event,gamma)*beta;
%compute the means at grid points in the considered event given estimates of b   
fu_event=g(u_event,gamma)*beta;
%compute the matrix of Euclidean distances among grid points and stations in 
%event under the consideration
c=sep(w_event,v_event);
%compute the covariance matrix among stations and the grid points given 
%estimates of tau2, sigma2 and h
switch cf
    case 'No'
Omega=TotalCovNon(theta,c);        
    case 'Exp'
Omega=TotalCovExp(theta,c);
    case 'SExp'
Omega=TotalCovSExp(theta,c); 
    case 'Matern1.5'
Omega=TotalCovMatern15(theta,c); 
end
%extract the covariance matrix among stations
n=size(x_event,1);
C=Omega(1:n,1:n);
%extract the covariance matrix between grid points and stations 
Sigma=Omega(n+1:end,1:n);
%extract the covariance matrix among grids
D=Omega(n+1:end,n+1:end);

mu=fu_event+Sigma/C*(y_event-fx_event);
V=D-Sigma/C*Sigma';
%Compute the predictions at grid points conditional on y_event
switch base
    case 'off'
z_hat=mu;
    case 'exp'
z_hat=exp(mu+0.5.*diag(V));
    case '10'
z_hat=exp(log(10).*mu+0.5*(log(10))^2.*diag(V));
end

end

function Omega=TotalCovNon(theta,c)
%This function computes the covariance matrix \Sigma without
%correlation function.  
Omega=theta(1)+theta(2).*eye(length(c));    
end

function Omega=TotalCovExp(theta,c)
%This function computes the covariance matrix \Sigma with exponential
%correlation function.  
Omega=theta(1)+theta(2).*exp(-c./theta(3));    
end

function Omega=TotalCovSExp(theta,c)
%This function computes the covariance matrix \Sigma with squared exponential
%correlation function.  
Omega=theta(1)+theta(2).*exp(-(c.^2)./(2*theta(3)^2));    
end

function Omega=TotalCovMatern15(theta,c)
%This function computes the covariance matrix \Sigma with Matern
%correlation function with v=1.5.
Omega=theta(1)+theta(2).*(1.+sqrt(3).*c./theta(3)).*exp(-sqrt(3).*c./theta(3));    
end

function c=sep(w_event,v_event)
%This function computes the speration distances between stations and grid
%points during an event.
position=[w_event;v_event];
m=size(position,1);
[K,L]=ndgrid(1:m);
K=K(:);
L=L(:);
pairs=[L,K];
X1=zeros(size(pairs,1),4);
for i=1:size(pairs,1)
X1(i,:)=[position(pairs(i,1),:), position(pairs(i,2),:)];
end
pairdis=deg2km(distance(X1(:,1),X1(:,2),X1(:,3),X1(:,4)));
c=vec2mat(pairdis, m);
end
