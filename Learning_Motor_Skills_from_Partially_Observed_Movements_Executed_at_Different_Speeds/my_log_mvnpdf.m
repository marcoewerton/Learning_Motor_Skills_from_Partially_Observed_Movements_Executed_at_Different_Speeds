function log_p = my_log_mvnpdf(x,mu,C)
    Sigma = chol(2*pi*C);
    logdetSigma_div2 = sum(log(diag(Sigma))); % logdetSigma divided by 2
    log_p =  -logdetSigma_div2 -(1/2)*(x-mu)*(C\(x-mu)');
end