function varargout = iv4(varargin)
    [lhs, rhs] = argn(0)
    
    plantData = varargin(1)
    orderData = varargin(2)
    na = orderData(1);nb = orderData(2)
    // arranging na ,nb,nk
    if size(orderData,"*") == 2 then
        nk = 1
    elseif size(orderData,'*') == 3 then
        nk = orderData(3)
    end
    nb1 = nb + nk - 1
    n = max(na, nb1)    
    // arranging the plant data
    if typeof(plantData) == 'constant' then
        Ts = 1;unitData = 'second'
    elseif typeof(plantData) == 'iddata' then
        Ts = plantData.Ts;unitData = plantData.TimeUnit
        plantData = [plantData.OutputData plantData.InputData]
    end
    noOfSample = size(plantData,'r')
    // finding the iv model
    ivTest = iv(plantData,[na nb nk]);
    // residual
    [aTemp,bTemp,cTemp] = pe(plantData,ivTest);
    Lhat = ar(aTemp,na+nb);
    x = sim(plantData(:,2),ivTest);
    yData = plantData(:,1);uData = plantData(:,2)
    Yf = filter(Lhat.a,Lhat.b,[plantData(:,1);zeros(n,1)]);
    phif = zeros(noOfSample,na+nb)
    psif = zeros(noOfSample,na+nb)
    // arranging samples of y matrix
    for ii = 1:na
        phif(ii+1:ii+noOfSample,ii) = -yData
        psif(ii+1:ii+noOfSample,ii) = -x
    end
    // arranging samples of u matrix
    for ii = 1:nb
        phif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
        psif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
    end
    // passing it through the filters
    for ii = 1:na+nb
        phif(:,ii) = filter(Lhat.a,Lhat.b,phif(:,ii));
        psif(:,ii) = filter(Lhat.a,Lhat.b,psif(:,ii));
    end
    lhs = psif'*phif
    lhsinv = pinv(lhs)
    theta = lhsinv * (psif)' * Yf
    ypred = (phif * theta)
    ypred = ypred(1:size(yData,'r'))
    e = yData - ypred
    sigma2 = norm(e)^2/(size(yData,'r') - na - nb)
    vcov = sigma2 * pinv((phif)' * phif)
    t = idpoly([1; theta(1:na)],[zeros(nk,1); theta(na+1:$)],1,1,1,Ts)
    
    // estimating the other parameters
    [temp1,temp2,temp3] = predict(z,t)
    [temp11,temp22,temp33] = pe(z,t)
    
    estData = calModelPara(temp1,temp11,na+nb)
    //pause
       t.Report.Fit.MSE = estData.MSE 
       t.Report.Fit.FPE = estData.FPE
    t.Report.Fit.FitPer = estData.FitPer
       t.Report.Fit.AIC = estData.AIC
      t.Report.Fit.AICc = estData.AICc
      t.Report.Fit.nAIC = estData.nAIC
       t.Report.Fit.BIC = estData.BIC
             t.TimeUnit = unitData
                    //sys = t
    varargout(1) = t
    //varargout(1) = idpoly([1; -theta(1:na)],[zeros(nk,1); theta(na+1:$)],1,1,1,Ts)
endfunction
