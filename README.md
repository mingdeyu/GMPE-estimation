# GMPE-estimation
This repository provides the MATLAB scripts used to estimate ground-motion prediction equations (GMPE) with spatial correlation. The scripts are produced as by-products of the forthcoming paper:

> Ming, D., Huang, C., Peters, G.W., and Galasso, C. An advanced estimation algorithm for ground-motion models with spatial correlation. *Bulletin of the Seismological Society of America* (forthcoming), 2018.

## General info
The current version of the estimation scripts only support the GMPE proposed by [Akkar and Bommer (2010)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/81/2/195/143661/empirical-equations-for-the-prediction-of-pga-pgv?redirectedFrom=fulltext):

<img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;f(\mathbf{b})=b_1&plus;b_2\,M_i&plus;b_3\,M_i^{2}&plus;(b_4&plus;b_5\,M_i)\log\sqrt{R_{ij}^2&plus;b_6^2}&plus;b_7\,S_{S,ij}&plus;b_8\,S_{A,ij}&plus;b_9\,F_{N,i}&plus;b_{10}\,F_{R,i}\,&space;\end{align*}" title="\small \begin{align*} f(\mathbf{b})=b_1+b_2\,M_i+b_3\,M_i^{2}+(b_4+b_5\,M_i)\log\sqrt{R_{ij}^2+b_6^2}+b_7\,S_{S,ij}+b_8\,S_{A,ij}+b_9\,F_{N,i}+b_{10}\,F_{R,i}\, \end{align*}" />

and four correlation functions:

* No spatial correlation:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=0&space;\end{align*}" title="\small \begin{align*} k(d)=0 \end{align*}" />

* Exponential:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=\exp\left(-\frac{d}{h}\right)&space;\end{align*}" title="\small \begin{align*} k(d)=\exp\left(-\frac{d}{h}\right) \end{align*}" />

* Matern with &nu;=1.5:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=\left(1&plus;\frac{\sqrt{3}d}{h}\right)\exp\left(-\frac{\sqrt{3}d}{h}\right)&space;\end{align*}" title="\small \begin{align*} k(d)=\left(1+\frac{\sqrt{3}d}{h}\right)\exp\left(-\frac{\sqrt{3}d}{h}\right) \end{align*}" />

* Squared exponential:

  <img src="https://latex.codecogs.com/svg.latex?\small&space;\begin{align*}&space;k(d)=\exp\left(-\frac{d^{2}}{2h^2}\right)&space;\end{align*}" title="\small \begin{align*} k(d)=\exp\left(-\frac{d^{2}}{2h^2}\right) \end{align*}" />

## Notes
The scripts can be easily extended to accommodate other forms of ground-motion prediction functions and correlation functions. If you would like to use the scripts on your own GMPE, please contact us and we are happy to help incorporate your GMPE into the scripts.

## Contents
The repository currently contains the following items:

* `scoring.m`: the MATLAB function that implements the Scoring estimation approach for GMPEs with spatial correlation;
* `data.mat`: the MATLAB data file that contains an example of ground-motion records as the inputs to ```scoring.m```;
* `Akkar2010.mat`: the MATLAB data file that contains the model parameter estimates of <b>b</b>, &tau; and &sigma; given by [Akkar and Bommer (2010)](https://pubs.geoscienceworld.org/ssa/srl/article-abstract/81/2/195/143661/empirical-equations-for-the-prediction-of-pga-pgv?redirectedFrom=fulltext); 
* `generator.m`: the MATLAB function that generates synthetic logarithmic intensity measures (IM).

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
%Set the initial value for the nonlinear coefficients (b<sub>6</sub> in this example)
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
[l,estimates,se,ci]=scoring(y,x,w,id,initial0,tol,cl,cf);

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

  - ```l```: the log-likelihood value at convergence;
  - ```estimates```: the estimates of all model parameters;
  - ```se```: the asymptotic standard error estimates of ML estimators of model parameters;
  - ```ci```: the confidence intervals (95% in this example) of all model parameters.

## A real data example
We will update this section soon...

## Built With
The scripts are built under MATLAB version R2018a.

## Contact
Please feel free to email me with any questions: 

Deyu Ming (deyu.ming.16@ucl.ac.uk).
