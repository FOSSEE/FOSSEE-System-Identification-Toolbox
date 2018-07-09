function varargout = iv(varargin)
    
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
    // arranging samples of y matrix
    for ii = 1:na
        phif(ii+1:ii+noOfSample,ii) = yData
        psif(ii+1:ii+noOfSample,ii) = x
    end
    // arranging samples of u matrix
    for ii = 1:nb
        phif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
        psif(ii+nk:ii+noOfSample+nk-1,ii+na) = uData
    end
    lhs = psif'*phif
    lhsinv = pinv(lhs)
    //pause
    theta = lhsinv * (psif)' * [yData;zeros(n,1)]
    ypred = (phif * theta)
    ypred = ypred(1:size(yData,'r'))
    e = yData - ypred
    sigma2 = norm(e)^2/(size(yData,'r') - na - nb)
    vcov = sigma2 * pinv((phif)' * phif)
    t = idpoly([1; -theta(1:na)],[zeros(nk,1); theta(na+1:$)],1,1,1,1)
    
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
             t.TimeUnit = unit
                    //sys = t
    varargout(1) = t
endfunction
