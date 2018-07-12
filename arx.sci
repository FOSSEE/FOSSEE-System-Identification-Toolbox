
// Estimates Discrete time ARX model
// A(q)y(t) = B(q)u(t) + e(t)
// Current version uses random initial guess
// 

// Authors: Ashutosh,Harpreet,Inderpreet
// Updated(12-6-16)

// Examples
//loadmatfile("data.mat")
//sys = arx(data,[2,2,1])
//sys  = 
//
//  A(z) = 1 - 1.3469229 z^-1 + 0.7420890 z^-2
//
//  B(z) = 1.3300766 z^-1 - 0.5726208 z^-2
//
//  Sampling Time = 1 seconds
//
//     MSE     FPE   FitPer        AIC      AICc    nAIC        BIC 
//  7.4091  7.4388  49.9726  9689.1801  194.2693  2.0067  9711.5838


function sys = arx(varargin)
	[lhs , rhs] = argn();	
	if ( rhs < 2 ) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "arx", rhs);
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
		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "arx");
		error(errmsg);
	end

	if (~isreal(z)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "arx");
		error(errmsg);
	end

	n = varargin(2)
	if (size(n,"*")<2| size(n,"*")>3) then
		errmsg = msprintf(gettext("%s: The order and delay matrix [na nb nk] should be of size [2 3]"), "arx");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [na nb nk] should be nonnegative integer number "), "arx");
		error(errmsg);
	end

	na = n(1); nb = n(2); //nk = n(3); //nf = n(4);
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
        e = YDATA - _objfunarx(UDATA,YDATA,p,na,nb,nk);
    endfunction
    tempSum = na+nb
    p0 = linspace(0.1,0.9,tempSum)';
    [var,errl] = lsqrsolve(p0,G,size(UDATA,"*"));
    err = (norm(errl)^2);
    opt_err = err;
	resid = G(var,[]);
    
    a = 1-poly([var(nb+1:nb+na)]',"q","coeff");
    b = poly([repmat(0,nk,1);var(1:nb)]',"q","coeff");
    a = (poly([1,-coeff(a)],'q','coeff'))
    t = idpoly(coeff(a),coeff(b),1,1,1,Ts)
    
    // estimating the other parameters
    [temp1,temp2,temp3] = predict(z,t)
    [temp11,temp22,temp33] = pe(z,t)
    
    estData = calModelPara(temp1,temp1,n(1)+n(2))
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
    
endfunction

function yhat = _objfunarx(UDATA,YDATA,x,na,nb,nk)
    x=x(:)
     q = poly(0,'q')
    tempSum = nb+na
    // making polynomials
    b = poly([repmat(0,nk,1);x(1:nb)]',"q","coeff");
    a = 1 - poly([x(nb+1:nb+na)]',"q","coeff")
    aSize = coeff(a);bSize = coeff(b)
    maxDelay = max([length(aSize) length(bSize)])
    yhat = [YDATA(1:maxDelay)]
    for k=maxDelay+1:size(UDATA,"*")
        tempB = 0
        for ii = 1:size(bSize,'*')
            tempB = tempB + bSize(ii)*UDATA(k-ii+1)
        end
        tempA = 0
        for ii = 1:size(aSize,"*")
            tempA = tempA + aSize(ii)*YDATA(k-ii)
        end
        yhat = [yhat; [ tempA+tempB ]];
    end
endfunction
