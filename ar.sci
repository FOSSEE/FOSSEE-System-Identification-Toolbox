
function sys = ar(varargin)
// Parameters Estimation of AR model using Input Output time-domain data
// 
// Calling Sequence
// sys = ar(ioData,[na])
// 
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// na : non-negative integer number specified as order of the polynomial A(z^-1)
// sys : idpoly type polynomial have estimated coefficients of A(z^-1) polynomials
// 
// Description
// Fit AR model on given input output data 
// The mathematical equation of the AR model 
// <latex>   
// begin{eqnarray}
// A(q)y(t) = e(t)
// end{eqnarray}
// </latex>
// It is SISO type model. It minimizes the sum of the squares of nonlinear functions using Levenberg-Marquardt algorithm.
// sys ,an idpoly type class, have different fields that contains estimated coefficients, sampling time, time unit and other estimated data in Report object.
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.5];b = [0 2 3];
//  model = idpoly(a,b,'Ts',0.1);
//  y = sim(u,model) + rand(length(u),1);
//  plantData = iddata(y,[],0.1);
//  sys = ar(plantData,[2])
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1]);
//  a = [1 0.5];b = [0 0.2 0.3];
//  model = idpoly(a,b,'Ts',0.1);
//  y = sim(u,model) + rand(length(u),1);
//  plantData = [y];
//  sys = ar(plantData,[2])
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar  


	[lhs , rhs] = argn();	
	if ( rhs < 2 ) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "ar", rhs);
			error(errmsg)
	end

	z = varargin(1)
    if typeof(z) == 'iddata' then
        Ts = z.Ts;unit = z.TimeUnit
        z = [z.OutputData z.InputData]
    elseif typeof(z) == 'constant' then
        Ts = 1;unit = 'seconds'
    end
	if ~iscolumn(z) then
		errmsg = msprintf(gettext("%s: time series output data only"), "ar");
		error(errmsg);
	end

	if (~isreal(z)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "ar");
		error(errmsg);
	end

	n = varargin(2)
	if (size(n,"*") ~=1 )then
		errmsg = msprintf(gettext("%s: order should be nonnegative integer number "), "ar");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [na] should be nonnegative integer number "), "ar");
		error(errmsg);
	end

	na = n; nb = 0; nk = 0; 
    //  storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = zeros(size(z,1),1)
    NDATA = size(UDATA,"*");
    function e = G(p,m)
        e = YDATA - _objfun(UDATA,YDATA,p,na,nb,nk);
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
    t = idpoly(coeff(a),1,1,1,1,Ts)
    
    //  estimating the other parameters
    [temp1,temp2,temp3] = predict([YDATA UDATA],t)
    [temp11,temp22,temp33] = pe([YDATA UDATA],t)
    
    estData = calModelPara(temp1,temp1,n(1))
    // pause
       t.Report.Fit.MSE = estData.MSE 
       t.Report.Fit.FPE = estData.FPE
    t.Report.Fit.FitPer = estData.FitPer
       t.Report.Fit.AIC = estData.AIC
      t.Report.Fit.AICc = estData.AICc
      t.Report.Fit.nAIC = estData.nAIC
       t.Report.Fit.BIC = estData.BIC
             t.TimeUnit = unit
                    sys = t
    // sys = idpoly(coeff(a),1,1,1,1,Ts)
//     sys.TimeUnit = unit
endfunction

function yhat = _objfun(UDATA,YDATA,x,na,nb,nk)
    x=x(:)
     q = poly(0,'q')
    tempSum = nb+na
    //  making polynomials
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
