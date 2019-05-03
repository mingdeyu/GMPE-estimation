%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%Copyright (C): Deyu MING
%Date: 2 May 2019
%Affiliation: Dept of Statistical Science at University College London
%Email: deyu.ming.16@ucl.ac.uk
%
% Reference: Ming, D., Huang, C., Peters, G.W., and Galasso, C. 
% An advanced estimation algorithm for ground-motion models with 
% spatial correlation. BSSA, 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=generator(x,w,id,beta,gamma,theta,g,n,cf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates n synthetic PGA data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% x: the covariates in the mean function f(X_ij);
%
% w: the covariates (i.e., the station coordinates) for the correlation function; 
%
% id: a vector containing values to differeniate earthquakes;
%
% beta: the chosen true values of linear coefficients in the mean function f(X_ij); 
%
% gamma: the chosen true values of nonlinear coefficients in the mean 
% function f(X_ij); 
%
% theta: the chosen true values of \tau2, \sigma2, and h;
%
% g: a function handle of the user-defined design matrix function for the
% linear coefficients;
%
% n: the number of the simulating trials to be generated;
%
% cf: the type of the covariance function: currently support 'exponential'
% type (Exp), Matern type with v=1.5 (Matern1.5) and 'squared exponential' type (SExp).
%
% OUTPUT:
% y: n synthetic PGA data with each column representing one of n trials;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=g(x,gamma)*beta;
c=sepdis(w,id);
switch cf
    case 'Exp'
Omega=CovExp(theta,c);
    case 'SExp'
Omega=CovSExp(theta,c);
    case 'Matern1.5'
Omega=CovMatern15(theta,c);
end
pre=repmat(f,1,n)';
y=mvnrnd(pre,Omega)';

end

function c=sepdis(position,id)
%This function computes the speration distances between stations for each
%earthquake.
%Inputs: id is the index labelling the event; position is the coordinates of
%stations
%Output: c is a cell containing the separation distance mateices for each
%event.
[g,event]=findgroups(id);
c=cell(1,length(event));
for j=1:length(event)
    subpos=position(g==j,:);
m=size(subpos,1);
[K,L]=ndgrid(1:m);
K=K(:);
L=L(:);
pairs=[L,K];
X1=zeros(size(pairs,1),4);
for i=1:size(pairs,1)
X1(i,:)=[subpos(pairs(i,1),:), subpos(pairs(i,2),:)];
end
pairdis=deg2km(distance(X1(:,1),X1(:,2),X1(:,3),X1(:,4)));
subdis=vec2mat(pairdis, m);
c{1,j}=subdis;
end
end

function Omega=CovExp(theta,c)
%This function computes the covariance matrix \Omega with exponential
%correlation function.  
n=length(c);
Omega_i=cell(1,n);
for i=1:n
Omega_i{1,i}=theta(1)+theta(2).*exp(-c{i}./theta(3));    
end
Omega=blkdiag(Omega_i{:});
end

function Omega=CovSExp(theta,c)
%This function computes the covariance matrix \Omega with squared exponential
%correlation function.  
n=length(c);
Omega_i=cell(1,n);
for i=1:n
Omega_i{1,i}=theta(1)+theta(2).*exp(-(c{i}.^2)./(2*theta(3)^2));    
end
Omega=blkdiag(Omega_i{:});
end

function Omega=CovMatern15(theta,c)
%This function computes the covariance matrix \Omega with Matern
%correlation function with v=1.5.
n=length(c);
Omega_i=cell(1,n);
for i=1:n
Omega_i{1,i}=theta(1)+theta(2).*(1.+sqrt(3).*c{i}./theta(3))...
    .*exp(-sqrt(3).*c{i}./theta(3));    
end
Omega=blkdiag(Omega_i{:});
end
