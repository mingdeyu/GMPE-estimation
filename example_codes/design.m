function B=design(x,gamma)
%This function gives the design matrix of \beta, the linear coefficients in
%the ground-motion prediction function of Akkar and Bommer (2010).
n=size(x,1);
B=[ones(n,1), x(:,1), x(:,1).^2, 0.5*log10(x(:,2).^2+gamma(1)^2),...
   0.5*x(:,1).*log10(x(:,2).^2+gamma(1)^2), x(:,3), x(:,4), x(:,5), x(:,6)];
end

%Change above to specify the design matrix of your own. When your GMPE
%has only linear coefficients, simply do not use gamma argument in coding
%the design matrix. However, keep the 'gamma' argument in the function definition
%on the first line.