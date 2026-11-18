function plotConvergence(res,titleText)
%PLOTCONVERGENCE relative residuals of the iterative solver
    if nargin < 2, titleText = 'Convergence'; end
    figure;
    semilogy(0:numel(res.resvec)-1,res.resvec,LineWidth=2);
    title(titleText,FontSize=20,Interpreter='latex');
    xlabel('iterations',FontSize=18,Interpreter='latex');
    ylabel('relative residual $\frac{|b-Ax_i|}{|b|}$',FontSize=18,Interpreter='latex');
end
