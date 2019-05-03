
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

function [varcomponents,semivar]=semivag(y,x,w,id,output,g,h0,deltah,cf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function draws the empirical semivarigram and the corresponding
% theorical semivariogram. It also estimate the range parameter h in 
% the correlation function.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% y: the logarithmic intensity measures;
%
% x: the covariates in the mean function f(X_ij);
% 
% w: the covariates (i.e., the station coordinates) in the correlation function; 
% 
% id: a vector containing the IDs of earthquakes; 
%
% output: the output structure produced from the scoring.m for the GMPE ignoring
% the spatial correlation;
% 
% g: a function handle of the user-defined design matrix function for the
% linear coefficients;
%
% h0: the initial value of h;
% 
% deltah: the seperating distance binwidth
% 
% cf: the type of the covariance function: currently support 'exponential'
% type (Exp), Matern type with v=1.5 (Matern1.5) and 'squared exponential'
% type (SExp); 
% 
% OUTPUTS:
% varcomponents: the ordinary nonlinear least squares estimates of intraevent 
% variance and range parameter;
% 
% semivar: the table of [d, semivariograms, N_d]:
%     1. d: a sequence of the seperating distances (the middle values of 
%        the bins);
%     2. semivariograms: the empirical semivariograms corresponding to d;
%     3. N_d: the number of pairs with seperating distances in the range of
%        the corresponding bins.
% Note: Following Jayaram and Baker (2009), the empirical semivariograms
% with d>100km are excluded. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=output.ParameterEstimates.beta;
gamma=output.ParameterEstimates.gamma;
varp0=[output.ParameterEstimates.theta(2);h0];
% Compute total residuals 
res=y-g(x,gamma)*beta;
% Group data based on earthquake event
[g, event]=findgroups(id);
% Count the station number for each event
n=zeros(size(event,1),1);
for i=1:size(event,1)
    n(i)=sum(g==i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This loop calculates the empirical semivariograms
for i=1:length(n)
    res_temp=res(g==i);
    w_temp=w(g==i,:);
    index=1;
    for j=1:n(i)
        for k=j+1:n(i)
            sep_dis(index)=deg2km(distance(w_temp(j,1),w_temp(j,2),w_temp(k,1),w_temp(k,2)));
            ressqr(index)=(res_temp(j)-res_temp(k)).^2;
            index=index+1;
        end
    end 
        gamma=table(repmat(i,length(ressqr),1),sep_dis',ressqr','VariableNames',{'Eventn_ID','InterStaDistance','ressqr'}); % res is intra+inter residual devided by intra J&Baker 2009
        if i==1
            gamma_out=gamma;
        else
            gamma_out=[gamma_out;gamma];
        end
        clear index sep_dis ressqr gamma
end
binnumber=round(100/(2*deltah));

for i=1:1:binnumber 
    index=find(gamma_out.InterStaDistance>(i-1)*deltah & gamma_out.InterStaDistance<=(i+1)*deltah);
    N_d(i)=size(index,1);
    semivariogram(i)=1/(2*N_d(i))*sum(gamma_out.ressqr(index));  % classical estimator Matheron, 1962
    clear index
end
d=(1:1:binnumber)'.*deltah.*2-deltah;
semivar=table(d,semivariogram',N_d','VariableNames',{'d_km','semivariogram','N_d'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part fits the parametric semivariogram model to the empirical semivariograms 
% by ordinary least squares
switch cf
    case 'Exp'
        spatailfcn=@(varp,d)varp(1).*(1-exp(-d./varp(2)));
    case 'SExp'
        spatailfcn=@(varp,d)varp(1).*(1-exp(-(d.^2)./(2.*varp(2).^2)));
    case 'Matern1.5'
        spatailfcn=@(varp,d)varp(1).*(1-(1+sqrt(3).*d./varp(2)).*exp(-sqrt(3).*d./varp(2)));
end

try
% fit semivariogram to find the intraevent variance, measurement variance/nugget effect
% and range parameter
    spatialmodelc=fitnlm(d,semivariogram,spatailfcn,varp0); 
% the estimates of intraevent variance, measurement variance/nugget effect
% and range parameter  
    varcomponents=spatialmodelc.Coefficients.Estimate; 
 
catch
disp('Enter grid search')
% in case the algorithm in function 'fitnlm' fails, the function
% 'gridsearch' is applied.
    
%Set initial interal lengths
interval_len=varp0.*10;
    varcomponents=gridsearch(d,semivariogram,interval_len,spatailfcn);
end
   
% plot empirical semivariogram
    figure;
    hold on;
    plot(d,semivariogram,'o');
    plot(d,spatailfcn(varcomponents,d));
    xlabel('d (km)');
    ylabel('Semivariogram');
    ax = gca;
    ax.FontSize = 14;
    switch cf
    case 'Exp'
        legend({'Empirical semivariogram','Exponential type'},'Location','southeast',...
            'FontSize',14)
    case 'SExp'
        legend({'Empirical semivariogram','Squared exponential type'},'Location',...
            'southeast','FontSize',14)
    case 'Matern1.5'
        legend({'Empirical semivariogram','Matérn type with \nu=1.5'},'Location',...
            'southeast','FontSize',14)
    end
    
end

function x_min=gridsearch(Xc,Yc,int_length,spatailfcn)
% This function find the range parameter by grid search technique
dim=length(int_length);
m1=50;
nr=20;
G1=linspace(0,int_length(1),m1);
G2=linspace(0,int_length(2),m1);
x_min=zeros(dim,1);
f_min=1e+10;

for r=1:nr
        for i=1:m1
        for j=1:m1
            f=sum((Yc-spatailfcn([G1(i);G2(j)],Xc)).^2);
            if f<f_min
                f_min=f;
                x_min=[G1(i);G2(j)];
            end
        end
        end
    int_length=int_length/2;
    G1=linspace(x_min(1)-int_length/2,x_min(1)+int_length/2,m1);
    G2=linspace(x_min(2)-int_length/2,x_min(2)+int_length/2,m1);
end
        
end
