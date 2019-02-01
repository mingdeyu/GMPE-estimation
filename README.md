<img src="https://raw.githubusercontent.com/mingdeyu/GMPE-estimation/master/logo/Logo.png" width="1000">

# GMPE-estimation

This repository provides the MATLAB scripts used to estimate ground-motion prediction equations (GMPE) with spatial correlation. The scripts are produced as by-products of the forthcoming paper:

> Ming, D., Huang, C., Peters, G.W., and Galasso, C. An advanced estimation algorithm for ground-motion models with spatial correlation. *Bulletin of the Seismological Society of America* (forthcoming), 2019.

## General info
The current version of the estimation scripts only support the GMPE proposed by [Akkar and Bommer (2010)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/81/2/195/143661/empirical-equations-for-the-prediction-of-pga-pgv?redirectedFrom=fulltext):

<img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;f(\mathbf{b})=b_1&plus;b_2\,M_i&plus;b_3\,M_i^{2}&plus;(b_4&plus;b_5\,M_i)\log\sqrt{R_{ij}^2&plus;b_6^2}&plus;b_7\,S_{S,ij}&plus;b_8\,S_{A,ij}&plus;b_9\,F_{N,i}&plus;b_{10}\,F_{R,i}\,&space;\end{align*}" title="\small \begin{align*} f(\mathbf{b})=b_1+b_2\,M_i+b_3\,M_i^{2}+(b_4+b_5\,M_i)\log\sqrt{R_{ij}^2+b_6^2}+b_7\,S_{S,ij}+b_8\,S_{A,ij}+b_9\,F_{N,i}+b_{10}\,F_{R,i}\, \end{align*}" />

and four correlation functions:

* No spatial correlation:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=0&space;\end{align*}" title="\small \begin{align*} k(d)=0 \end{align*}" />

* Exponential:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=\exp\left(-\frac{d}{h}\right)&space;\end{align*}" title="\small \begin{align*} k(d)=\exp\left(-\frac{d}{h}\right) \end{align*}" />

* Mat&eacute;rn with &nu;=1.5:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=\left(1&plus;\frac{\sqrt{3}d}{h}\right)\exp\left(-\frac{\sqrt{3}d}{h}\right)&space;\end{align*}" title="\small \begin{align*} k(d)=\left(1+\frac{\sqrt{3}d}{h}\right)\exp\left(-\frac{\sqrt{3}d}{h}\right) \end{align*}" />

* Squared exponential:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=\exp\left(-\frac{d^{2}}{2h^2}\right)&space;\end{align*}" title="\small \begin{align*} k(d)=\exp\left(-\frac{d^{2}}{2h^2}\right) \end{align*}" />

## Notes
The scripts can be easily extended to accommodate other forms of ground-motion prediction functions and correlation functions. If you would like to use the scripts on your own GMPE, please contact us and we are happy to help incorporate your GMPE into the scripts.

## Contents
The repository currently contains the following items:

* `Akkar2010.mat`: the MATLAB data file that contains the model parameter estimates of <b>b</b>, &tau; and &sigma; given by [Akkar and Bommer (2010)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/81/2/195/143661/empirical-equations-for-the-prediction-of-pga-pgv?redirectedFrom=fulltext);
* `data.mat`: the MATLAB data file that contains an example of ground-motion records as the inputs to ```scoring.m```; 
* `eventfinder.m`: the MATLAB function that extracts information of a specified event from the ground-motion database;
* `generator.m`: the MATLAB function that generates synthetic logarithmic intensity measures (IM).
* `grid_data.mat`: the MATLAB data file that contains the covariate information at prediction locations;
* `Italyborder.mat`: the MATLAB data file that contains the coordinates of the border of Italy;
* `prediction.m`: the MATLAB function that make linear predictions on the specified locations given the observations at recording sites;
* `scoring.m`: the MATLAB function that implements the Scoring estimation approach for GMPEs with spatial correlation;
* `semivag.m`: the MATLAB function that plots the empirical semivariogram together with the fitted theoretical semivariogram models;
* `shakemap.m`: the MATLAB function that draws the shake map based on the produced predictions.

## A synthetic example
This example illustrates a step-by-step instruction on how to use the Scoring estimation approach to estimate a GMPE with spatial correlation based on a synthetic dataset of logarithmic PGA. 

* **Step 1 - Data preparation** 

  In this step, we load and assign data required for the estimation.

```
%load data
load('data.mat')

%Assign data to covariates in the ground-motion prediction function given by Akkar and Bommer (2010).
x=[data.Mw,data.JB_dist,data.Ss,data.Sa,data.Fn,data.Fr];

%Assign data to covariates (lat and lon) required in the correlation function
w=[data.st_latitude,data.st_longitude];

%Set earthquake ID
id=data.event_id;
```

* **Step 2 - Synthetic data generation** 

  To generate synthetic logarithmic IMs (PGA in this case), we need to decide the type of correlation function to, choose values of model parameters, and pick the number of synthetic logarithmic PGA dataset to be generated.

``` 
%Choose correlation function type ('No','Exp' or 'Matern1.5' or 'SExp')
cf='Exp';

%Assign parameter values given by Akkar and Bommer (2010)
%para(1:10) correspond to b, para(11) and para(12) correspond to tau and sigma, respectively
load('Akkar2010.mat')
para=Akkar2010; 

%Assign a value to the range parameter h in the exponential correlation function
h=11.5; 

%Generate one set of log(IM) (i.e., log(PGA))
n=1;
y=generator(x,w,id,para,h,n,cf);
```

* **Step 3 - Estimation** 

  The estimation of GMPE is then can be implemented by executing the following scripts:

```
%Set the initial value for the nonlinear coefficients (b_6 in this example)
gamma0=para(6)-1;

%Set initial values for the inter- and intraevent variances and the range parameter h
theta0=[para(11)^2-0.001,para(12)^2-0.01,h-2]';
initial0=[gamma0;theta0];

%In this example, we perturb the parameter values used to generate synthetic data as the initial values for the estimation

%Set the tolerance level for the optimization algorithm
tol=0.0001;

%Set the confidence level (95% in this example) to be used to construct confidence intervals for the parameter estimators
cl=95;

%Begin the estimation
[info,l,estimates,se,ci]=scoring(y,x,w,id,initial0,tol,cl,cf);

%%%%%%%%%%%%COMMAND WINDOW OUTPUTS%%%%%%%%%%%%
  Iteration 0 Loglikelihood   107.7564
  Iteration 1 Loglikelihood   120.1908
  Iteration 2 Loglikelihood   120.3800
  Iteration 3 Loglikelihood   120.3845
  Iteration 4 Loglikelihood   120.3847
  Converged!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Note: your outputs will differ from the above depending on your synthetic data. 

```

* **Step 4 - Extraction of results** 

  The estimation results can be extracted from the following four variables:

  - ```info```: the AIC and BIC values;
  - ```l```: the log-likelihood value at convergence;
  - ```estimates```: the estimates of all model parameters;
  - ```se```: the asymptotic standard error estimates of ML estimators of model parameters;
  - ```ci```: the confidence intervals (95% in this example) of all model parameters.

## A real data example
In this section, we demonstrate in a real data example on how to use the Scoring estimation approach to estimate the GMPE proposed by Akkar and Bommer (2010) with spatial correlation. 

* **Step 1 - Initial estimation** 

  In this step, we estimate the GMPE proposed by Akkar and Bommer (2010) without considering the spatial correlation. The results of this initial estimation will provide information to preliminarily estimate spatial correlation function in Step 2.

```
%load data
load('data.mat')

%Assign data to covariates in the ground-motion prediction function given by Akkar and Bommer (2010).
x=[data.Mw,data.JB_dist,data.Ss,data.Sa,data.Fn,data.Fr];

%Assign data to covariates (lat and lon) required in the correlation function
w=[data.st_latitude,data.st_longitude];

%Set earthquake ID
id=data.event_id;

%Log transform the PGAs and assign them to the dependent variable
y=log10(data.rotD50_pga);

%Choose correlation function type ('No','Exp' or 'Matern1.5' or 'SExp')
cf='No';

%Assign parameter values given by Akkar and Bommer (2010)
%para(1:10) correspond to b, para(11) and para(12) correspond to tau and sigma, respectively
load('Akkar2010.mat')
para=Akkar2010; 

%Set the initial value for the nonlinear coefficients (b_6 in this example) as given by Akkar and Bommer (2010)
gamma0=para(6);

%Set initial values for the inter- and intraevent variances as given by Akkar and Bommer (2010)
theta0=[para(11)^2,para(12)^2]';
initial0=[gamma0;theta0];

%Set the tolerance level for the optimization algorithm
tol=0.0001;

%Set the confidence level (95% in this example) to be used to construct confidence intervals for the parameter estimators
cl=95;

%Begin the estimation
[info,l,estimates,se,ci]=scoring(y,x,w,id,initial0,tol,cl,cf);

%%%%%%%%%%%%COMMAND WINDOW OUTPUTS%%%%%%%%%%%%
  Iteration 0 Loglikelihood  -877.7142
  Iteration 1 Loglikelihood  -728.0676
  Iteration 2 Loglikelihood  -711.1419
  Iteration 3 Loglikelihood  -709.8849
  Iteration 4 Loglikelihood  -709.8144
  Iteration 5 Loglikelihood  -709.8132
  Iteration 6 Loglikelihood  -709.8132
  Converged!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
``` 

* **Step 2 - Preliminary estimation of the range paramter** 

  In this step, we use the estimates obtained in Step 1 to obtain the preliminary estimate of the range parameter via the semivariogram.

```
%Set an initial value for the range parameter h
h0=40;

%Set the binwidth for the empirical semivariogram
deltah=0.5;

%Plot empirical semivariogram with the fitted exponential semivariograms
[varcomponents_exp,semivar]=semivag(y,x,w,id,estimates,h0,deltah,'Exp');
```
  This will produce the following figure:

 <img src="https://raw.githubusercontent.com/mingdeyu/GMPE-estimation/master/tutorial_figures/exp.png" width="500">

* **Step 3 - Estimation** 

  In this step, we use the preliminary estimate of h obtained in Step 2 to implement the one-stage estimation algorithm.

```
%Set initial value for the nonlinear coefficient b_6 as the one estimated from Step 1 
gamma0=estimates(6);

%Set initial value for the interevent variance tau^2 as the one estimated from Step 1
tau20=estimates(11);

%Set initial values for the intraevent variance sigma^2 and the range parameter h as the ones estimated from Step 2
sigma20_exp=varcomponents_exp(1);
h0_exp=varcomponents_exp(2);

%Combine initial values into one vector
initial0_exp=[gamma0;tau20;sigma20_exp;h0_exp];

%Begin estimation
[info_exp,ll_exp,estimates_exp,se_exp,ci_exp]=scoring(y,x,w,id,initial0_exp,tol,cl,'Exp');

%%%%%%%%%%%%COMMAND WINDOW OUTPUTS%%%%%%%%%%%%
      .             .             .
      .             .             .
      .             .             .
  Iteration 32 Loglikelihood  -675.1151
  Iteration 33 Loglikelihood  -675.1149
  Iteration 34 Loglikelihood  -675.1147
  Iteration 35 Loglikelihood  -675.1146
  Iteration 36 Loglikelihood  -675.1144
  Iteration 37 Loglikelihood  -675.1143
  Converged!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

  Step 1 to 3 can be repeated for other choices of correlation functions, e.g., Mat&eacute;rn with &nu;=1.5 and squared exponential. However, the preliminary estimate of the range parameter obtained in Step 2 can be very biased as demonstrated in the paper, giving an unreasonable initial value to h, which in turn will cause ill-conditioning problems (for the squared exponential type, particularly) to the estimation algorithm. In such case, users are suggested to try a smaller initial value for h and re-run the estimation algorithm.  

  Though we will not go through in this example, users of the scripts are suggested to always check the statistical significances of model parameters via the confidence intervals produced by the estimation algorithm. This will help users decide whether some varaibles need to be removed, whether the model should be modified, whether more data are required to improve the statistical significance, etc. For example, in this real data case the magnitude and fault types are not statistical significant. This is possibly because the information of the magnitude and fault types are given per event rather than per record, giving insufficient evidences to show they are statistically siginificant. Therefore, more events should be included in the training data. However, to prevent from the curse of dimensionality, it might be a good idea to include small events which have small numbers of records to the data. 

* **Step 4 - Prediction** 

  In this step, we show how to draw a ShakeMap for a given event based on the GMPE estimated in Step 3.
```
%Define the event under consideration
eventid='IT-1997-0137';

%Extract information of the event from the ground-motion database
[y_event,x_event,w_event]=eventfinder(eventid,id,y,x,w);

%load information at grid points (i.e., prediction locations)
load('grid_data.mat')

%Assign data to covariates in the ground-motion prediction function given by Akkar and Bommer (2010).
u_event=[grid_data.Mw,grid_data.JB_dist,grid_data.Ss,grid_data.Sa,grid_data.Fn,grid_data.Fr];

%Assign data to covariates (lat and lon) required in the correlation function
v_event=[grid_data.st_latitude,grid_data.st_longitude];

%Make predictions on grid points
z_hat_exp=prediction(y_event,x_event,u_event,w_event,v_event,estimates_exp,'Exp','off');

%load event map info
load('Italyborder.mat')

%Define the region of the shakemap
latRange=[40.68,45.18];
lonRang=[10.68,15.18];

%Draw the shakemap without the log-scaled color (0) bar
lim=[min(z_hat_exp),max(z_hat_exp)]
f=shakemap(latRange,lonRang,Italyborder,w_event,v_event,z_hat_exp,0,lim);
title('ShakeMap of M_W 5.6','fontsize',10);
hcb=colorbar;
title(hcb,'log_{10}(PGA)','fontsize',6);
```

 The produced shake map in terms of log<sub>10</sub>(PGA) is shown below. 

<img src="https://raw.githubusercontent.com/mingdeyu/GMPE-estimation/master/tutorial_figures/shakemap_logpga.png" width="500">

 One can also plot the shake map in terms of PGA by:

```
%Make predictions on grid points and set the base to '10'
z_hat_exp_delog=prediction(y_event,x_event,u_event,w_event,v_event,estimates_exp,'Exp','10');

%Draw the shakemap with log-scaled (1) color bar
lim=[min(z_hat_exp_delog),max(z_hat_exp_delog)]
f=shakemap(latRange,lonRang,Italyborder,w_event,v_event,z_hat_exp_delog,1,lim);
title('ShakeMap of M_W 5.6','fontsize',10);
hcb=colorbar;
title(hcb,'PGA (cm/s^2)','fontsize',6);
```

The resulting shake map is shown below:

<img src="https://raw.githubusercontent.com/mingdeyu/GMPE-estimation/master/tutorial_figures/shakemap_pga.png" width="500">

## Built with
The scripts are built under MATLAB version R2018a.

## Contact
Please feel free to email me with any questions: 

Deyu Ming (deyu.ming.16@ucl.ac.uk).
