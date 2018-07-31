function varargout = iv4(varargin)
// Parameters Estimation of IV4 model by four stage instrumental variable method
// 
// Calling Sequence
// sys = iv(ioData,[na nb nk])
// 
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// na : non-negative integer number specified as order of the polynomial A(z^-1)
// nb : non-negative integer number specified as order of the polynomial B(z^-1)+1
// nk : non-negative integer number specified as input output delay, Default value is 1
// sys : idpoly type polynomial have estimated coefficients of A(z^-1) and B(z^-1) polynomials
// 
// Description
// Fit IV4 model on given input output data 
// The structure of sys is ARX type.The mathematical equation is given here 
// <latex>   
// begin{eqnarray}
// A(q)y(n) = B(q)u(n-k) + e(t)
// end{eqnarray}
// </latex>
// IV4 model is SISO type model. It is unaffected by color of the noise. Four steps used in IV4 model design. First step is the generation of the ARX model.   
// Second step uses the ARX model to generate the instrument variable matrix.Next steps uses the residual to generate a higher order model coefficient. 
// In final step uses the AR model coefficient to filter the input and output data and feed it to the IV model. 
// sys ,an idpoly type class, have different fields that contains estimated coefficients, sampling time, time unit and other estimated data in Report object.
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// ioData = iddata(y,u,0.1)
// sys = iv4(ioData,[2,2,1])
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// ioData = [y,u]
// sys = iv4(ioData,[2,2,1])
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar  

    [lhs, rhs] = argn(0)
    plantData = varargin(1)
    orderData = varargin(2)
    na = orderData(1);nb = orderData(2)
    //  arranging na ,nb,nk
    if size(orderData,"*") == 2 then
        nk = 1
    elseif size(orderData,'*') == 3 then
        nk = orderData(3)
    end
    nb1 = nb + nk - 1
    n = max(na, nb1)    
    //  arranging the plant data
    if typeof(plantData) == 'constant' then
        Ts = 1;unitData = 'second'
    elseif typeof(plantData) == 'iddata' then
        Ts = plantData.Ts;unitData = plantData.TimeUnit
        plantData = [plantData.OutputData plantData.InputData]
    end
    noOfSample = size(plantData,'r')
    //  finding the iv model
    ivTest = iv(plantData,[na nb nk]);
    //  residual
    [aTemp,bTemp,cTemp] = pe(plantData,ivTest);
    Lhat = ar(aTemp,na+nb);
    x = sim(plantData(:,2),ivTest);
    yData = plantData(:,1);uData = plantData(:,2)
    Yf = filter(Lhat.a,Lhat.b,[plantData(:,1);zeros(n,1)]);
    phif = zeros(noOfSample,na+nb)
    psif = zeros(noOfSample,na+nb)
    //  arranging samples of y matrix
    for ii = 1:na
        phif(ii+1:ii+noOfSample,ii) = -yData
        psif(ii+1:ii+noOfSample,ii) = -x
    end
    //  arranging samples of u matrix
    for ii = 1:nb
        phif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
        psif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
    end
    //  passing it through the filters
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
    
    //  estimating the other parameters
    [temp1,temp2,temp3] = predict(plantData,t)
    [temp11,temp22,temp33] = pe(plantData,t)
    
    estData = calModelPara(temp1,temp11,na+nb)
    // pause
       t.Report.Fit.MSE = estData.MSE 
       t.Report.Fit.FPE = estData.FPE
    t.Report.Fit.FitPer = estData.FitPer
       t.Report.Fit.AIC = estData.AIC
      t.Report.Fit.AICc = estData.AICc
      t.Report.Fit.nAIC = estData.nAIC
       t.Report.Fit.BIC = estData.BIC
             t.TimeUnit = unitData
                    // sys = t
    varargout(1) = t
    // varargout(1) = idpoly([1; -theta(1:na)],[zeros(nk,1); theta(na+1:$)],1,1,1,Ts)
endfunction
