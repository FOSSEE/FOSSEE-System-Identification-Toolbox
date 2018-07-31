function varargout = iv(varargin)
    
// Parameters Estimation of IV model by arbitrary instrumental variable method
// 
// Calling Sequence
// sys = iv(ioData,[na nb nk])
// sys = iv(ioData,[na nb nk],instData)
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// na : non-negative integer number specified as order of the polynomial A(z^-1)
// nb : non-negative integer number specified as order of the polynomial B(z^-1)+1
// nk : non-negative integer number specified as input output delay, Default value is 1
// instData : arbitrary instrument variable matrix. The size of instriment variable must be equal to the size of outputData
// sys : idpoly type polynomial have estimated coefficients of A(z^-1) and B(z^-1) polynomials
// 
// Description
// Fit IV model on given input output data 
// The mathematical equation of the IV model 
// <latex>   
// begin{eqnarray}
// A(q)y(n) = B(q)u(n-k) + e(t)
// end{eqnarray}
// </latex>
// It is SISO type model. Instrument variable method is use to estimate the cofficients. If user does not provide the arbitrary instrument variable matrix
// then the program extracte it by using ARX method.
// 
// sys ,an idpoly type class, have different fields that contains estimated coefficients, sampling time, time unit and other estimated data in Report object.
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// ioData = iddata(y,u,0.1)
// sys = arx(ioData,[2,2,1])
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// ioData = [y,u]
// sys = iv(ioData,[2,2,1])
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar
  
	[lhs , rhs] = argn(0);	
	if ( rhs < 2 || rhs > 3) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2 or 3"), "iv", rhs);
			error(errmsg)
	end
    plantData = varargin(1)
    if typeof(plantData) == 'iddata' then
        Ts = plantData.Ts;unit = plantData.TimeUnit
        plantData = [plantData.OutputData plantData.InputData]
    elseif typeof(plantData) == 'constant' then
        Ts = 1;unit = 'seconds'
    end
	if ((~size(plantData,2)==2) & (~size(plantData,1)==2)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "iv");
		error(errmsg);
	end

	if (~isreal(plantData)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "arx");
		error(errmsg);
	end

	n = varargin(2)
	if (size(n,"*")<2| size(n,"*")>3) then
		errmsg = msprintf(gettext("%s: The order and delay matrix [na nb nk] should be of size [2 or 3]"), "iv");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [na nb nk] should be nonnegative integer number "), "iv");
		error(errmsg);
	end
    na = n(1);nb = n(2)
    if size(n,'*') == 2 then
        nk = 1
    elseif size(n,'*') == 3 then
        nk = n(3)
    end
    yData = plantData(:,1)
    uData = plantData(:,2)
    noOfSample = size(plantData,'r')
    nb1 = nb + nk - 1
    n = max(na, nb1)    
    if rhs == 3 then
        if typeof(varargin(3)) <> 'constant' then
            errmsg = msprintf(gettext("%s: Incompatible last input argument. "), "iv");
            error(errmsg)
        elseif size(varargin(3),'r') <> size(uData,'r') then
            errmsg = msprintf(gettext("%s: number of samples of output must be equal to the dimensions of plant data "), "iv");
            error(errmsg);
        end
        x = varargin(3)
    elseif rhs == 2
        arxModel = arx(plantData,[na nb nk])
        x = sim(uData,arxModel)
    end
    phif = zeros(noOfSample,na+nb)
    psif = zeros(noOfSample,na+nb)
    //  arranging samples of y matrix
    for ii = 1:na
        phif(ii+1:ii+noOfSample,ii) = yData
        psif(ii+1:ii+noOfSample,ii) = x
    end
    //  arranging samples of u matrix
    for ii = 1:nb
        phif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
        psif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
    end
    lhs = psif'*phif
    lhsinv = pinv(lhs)
    // pause
    theta = lhsinv * (psif)' * [yData;zeros(n,1)]
    ypred = (phif * theta)
    ypred = ypred(1:size(yData,'r'))
    e = yData - ypred
    sigma2 = norm(e)^2/(size(yData,'r') - na - nb)
    vcov = sigma2 * pinv((phif)' * phif)
    t = idpoly([1; -theta(1:na)],[zeros(nk,1); theta(na+1:$)],1,1,1,1)
    
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
             t.TimeUnit = unit
                    // sys = t
    varargout(1) = t
endfunction
