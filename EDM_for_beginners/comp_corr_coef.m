% To test for convergence in convergent cross mapping (CCM), two approaches 
% are widely used. First, the convergence can be tested by investigating how 
% the cross- mapping skill changes with respect to the library size (e.g., 
% trend or increment). For example, one can consider the following two 
% statistical criteria: (1) testing the existence of a significant monotonic 
% increasing trend in q(L) using Kendall?s s test, and (2) testing the 
% significance of the improvement in q(L) by Fisher?s Dq Z test, which checks
% whether the cross-mapping skill obtained under the maximal library length 
% (q(Lmax)) is significantly higher than that obtained using the minimal library 
% length (q(L0)). The convergence of CCM is deemed significant when both 
% Kendall?s s test and Fisher?s Dq Z test are significant. 

% This function compares if two correlation coefficients are significantly 
% different. 
% The correlation coefficients were tansfered to z scores using fisher's r 
% to z transformation. 
% ref: http://core.ecu.edu/psyc/wuenschk/docs30/CompareCorrCoeff.pdf 
%-------------------------------------------------------------------------- 
% Inputs: (1) r1: correlation coefficient of the first correlation (2) r2: 
% correlation coefficient of the second correlation (3) n1: number of 
% samples used to compute the first correlation (4) n2: number of samples 
% used to compute the second correlation 
%-------------------------------------------------------------------------- 
% Output: (1) p: p value, the probability that H0 (the correlation 
% coefficiets are not different) is correct 
%-------------------------------------------------------------------------- 


r1=; 
r2=;
n1=;
n2=;
p = compare_correlation_coefficients(r1,r2,n1,n2); 


