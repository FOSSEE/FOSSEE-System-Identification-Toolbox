function varargout = fitValue(varargin)
    yData = varargin(1)
    ySys = varargin(2)
    fitValueNum = 0
    fitValueDen = 0
    sampleLength = size(yData,'r')
    yDataMean = mean(yData)
    for ii = 1:sampleLength
        fitValueNum = fitValueNum + (yData(ii)-ySys(ii))^2
        fitValueDen = fitValueDen + (yData(ii)-yDataMean)^2
    end
    fit = 100*(1-sqrt(fitValueNum/fitValueDen))
    varargout(1) = fit
endfunction
