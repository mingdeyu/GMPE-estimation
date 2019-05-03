function A=Bg(x,gamma)
% This function computes the grident of design matrix B wrt nonlinear coefficient
% gamma, B_{\gamma}
A=cell(1,length(gamma));

n=9; %Replace here by the number of linear coefficients in your GMPE!

% The following block creates the grident of design matrix B wrt the nonlinear 
% coefficient b_6 in the ground-motion prediction fct of Akkar&Bommer
% (2010)
%%%%%%%%%%%%%%%%%%%%%%%
%       Block 1 
%%%%%%%%%%%%%%%%%%%%%%%
M=zeros(size(x,1),n);
M(:,4)=gamma/log(10)./(x(:,2).^2+gamma^2);
M(:,5)=gamma/log(10).*x(:,1)./(x(:,2).^2+gamma^2);
A{1}=M;
%%%%%%%%%%%%%%%%%%%%%%%

%Since in Akkar and Bommer (2010) there is only one nonlinear coefficient,
%only one block (see above) is needed. However, if there are more than one 
%nonlinear coefficients, then the users need to create more blocks similar
%to the above, where each block computes the grident of design matrix B wrt
%each nonlinear coefficient in the gamma (i.e., each nonlinear coefficient
%in the ground-motion prediction function, e.g.,

%%%%%%%%%%%%%%%%%%%%%%%
%       Block 2 
%%%%%%%%%%%%%%%%%%%%%%%
% M=zeros(size(x));
% TYPE YOUR CODES HERE
% A{2}=M;

%%%%%%%%%%%%%%%%%%%%%%%
%       Block 3 
%%%%%%%%%%%%%%%%%%%%%%%
% M=zeros(size(x));
% TYPE YOUR CODES HERE
% A{3}=M;

% More blocks go under here......

end