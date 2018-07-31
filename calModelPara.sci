function varargout = calModelPara(varargin)
    y = varargin(1)+varargin(2)
    N = size(y,'r')
    ek = varargin(2)
    np = varargin(3)
    //pause
    mse = sum(ek.^2)/N;//disp(mse)
    fpe = mse * (1 + np/N)/(1 - np/N);//disp(fpe)
    nrmse = 1 - sqrt(sum(ek^2))/sqrt(sum((y - mean(y))^2));//disp(nrmse)
    AIC = N * log(mse) + 2 * np + N * size(y,'c') * (log(2 * %pi) + 1);//disp(AIC)
    AICc = AIC * 2 * np * (np + 1)/(N - np - 1);//disp(AICc)
    nAIC = log(mse) + 2 * np/N;//disp(nAIC)
    BIC = N * log(mse) + N * size(y,'c') * (log(2 * %pi) + 1) + np * log(N);//disp(BIC)
   // pause
    varargout(1) = struct('MSE',mse,'FPE',fpe,'FitPer',nrmse*100,'AIC',AIC,'AICc',AICc,'nAIC',nAIC,'BIC',BIC)
endfunction
