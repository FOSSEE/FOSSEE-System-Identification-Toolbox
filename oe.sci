// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Authors: Harpreet, Ashutosh
// Organization: FOSSEE, IIT Bombay

function sys =  oe(varargin)

// Estimates Discrete time BJ model
// y(t) = [B(q)/F(q)]u(t) + [C(q)/D(q)]e(t)
// Current version uses random initial guess
// Need to get appropriate guess from OE and noise models

	[lhs , rhs] = argn();	
	if ( rhs < 2 ) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "oe_2", rhs);
			error(errmsg)
	end


	z = varargin(1)
    if typeof(z) == 'iddata' then
        Ts = z.Ts;unit = z.TimeUnit
        z = [z.OutputData z.InputData]
    elseif typeof(z) == 'constant' then
        Ts = 1;unit = 'seconds'
    end
	if ((~size(z,2)==2) & (~size(z,1)==2)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "oe_2");
		error(errmsg);
	end

	if (~isreal(z)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "oe_2");
		error(errmsg);
	end
//
	n = varargin(2)
	if (size(n,"*")<2| size(n,"*")>3) then
		errmsg = msprintf(gettext("%s: The order and delay matrix [nb nf nk] should be of size [2 4]"), "oe_2");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [nb nf nk] should be nonnegative integer number "), "oe_2");
		error(errmsg);
	end
//
	nb= n(1); nf = n(2);
//	
	if (size(n,"*") == 2) then
		nk = 1
	else
		nk = n(3);
	end

    // storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = z(:,2);
    NDATA = size(UDATA,"*");
    function e = G(p,m)
        e = YDATA - _objoefun(UDATA,YDATA,p,nf,nb,nk);
    endfunction
    tempSum = nf+nb
    p0 = linspace(0.04,0.041,tempSum)';
    [var,errl] = lsqrsolve(p0,G,size(UDATA,"*"));
    err = (norm(errl)^2);
    opt_err = err;
	resid = G(var,[]);
    f = poly([1; var(nb+1:nb+nf)],"q","coeff");
    b = poly([repmat(0,nk,1);var(1:nb)]',"q","coeff");
    t = idpoly(1,coeff(b),1,1,coeff(f),Ts)
    
    // estimating the other parameters
    [temp1,temp2,temp3] = predict(z,t)
    [temp11,temp22,temp33] = pe(z,t)
    
    estData = calModelPara(temp1,temp11,n(1)+n(2))
    //pause
       t.Report.Fit.MSE = estData.MSE 
       t.Report.Fit.FPE = estData.FPE
    t.Report.Fit.FitPer = estData.FitPer
       t.Report.Fit.AIC = estData.AIC
      t.Report.Fit.AICc = estData.AICc
      t.Report.Fit.nAIC = estData.nAIC
       t.Report.Fit.BIC = estData.BIC
             t.TimeUnit = unit
                    sys = t
    //sys = t
    //sys.TimeUnit = unit
endfunction

function yhat = _objoefun(UDATA,YDATA,x,nf,nb,nk)
    x=x(:)
     q = poly(0,'q')
    tempSum = nb+nf
    // making polynomials
    b = poly([repmat(0,nk,1);x(1:nb)]',"q","coeff");
    f =  -1*poly([x(nb+1:nb+nf)]',"q","coeff")
    fSize = coeff(f);bSize = coeff(b)
    maxDelay = max([length(fSize) length(bSize)])
    yhat = [YDATA(1:maxDelay)]
    for k=maxDelay+1:size(UDATA,"*")
        tempB = 0
        for ii = 1:size(bSize,'*')
            tempB = tempB + bSize(ii)*UDATA(k-ii+1)
        end
        tempF = 0
        for ii = 1:size(fSize,"*")
            tempF = tempF + fSize(ii)*yhat(k-ii)
        end
        yhat = [yhat; [ tempB + tempF ]];
    end
endfunction
