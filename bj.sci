
function sys = bj(varargin)
    
// Parameters Estimation of BJ(Box-Jenkins) model using Input Output time-domain data
// 
// Calling Sequence
// sys = bj(ioData,[nb nc nd nf nk])
// 
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// nb : non-negative integer number specified as order of the polynomial B(z^-1)+1
// nc : non-negative integer number specified as order of the polynomial C(z^-1)
// nd : non-negative integer number specified as order of the polynomial D(z^-1)
// nf : non-negative integer number specified as order of the polynomial f(z^-1)
// nk : non-negative integer number specified as input output delay, Default value is 1
// sys : idpoly type polynomial have estimated coefficients of B(z^-1),C(z^-1),D(z^-1) and f(z^-1) polynomials
// 
// Description
// Fit BJ model on given input output data 
// The mathematical equation of the BJ model 
// <latex>   
// begin{eqnarray}
// y(n) = \frac {B(q)}{D(q)}u(n) + \frac {C(q)}{D(q)}e(t)
// end{eqnarray}
// </latex>
// It is SISO type model. It minimizes the sum of the squares of nonlinear functions using Levenberg-Marquardt algorithm.
// sys ,an idpoly type class, have different fields that contains estimated coefficients, sampling time, time unit and other estimated data in Report object.
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.5];b = [0 2 3];
//  model = idpoly(a,b,'Ts',0.1)
//  y = sim(u,model) + rand(length(u),1)
//  ioData = iddata(y,u,0.1)
//  sys = bj(ioData,[2,2,2,2,1])
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.5];b = [0 2 3];
//  model = idpoly(a,b,'Ts',0.1)
//  y = sim(u,model) + rand(length(u),1)
//  ioData = [y,u]
//  sys = bj(ioData,[2,2,2,2,1])
// 
// Authors
// Ashutosh Kumar Bhargava, Harpreet,Inderpreet  

	[lhs , rhs] = argn();	
	if ( rhs < 2 ) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "bj", rhs);
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
		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "bj");
		error(errmsg);
	end

	if (~isreal(z)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "bj");
		error(errmsg);
	end

	n = varargin(2)
	if (size(n,"*")<4| size(n,"*")>5) then
		errmsg = msprintf(gettext("%s: The order and delay matrix [nb nc nd nf nk] should be of size [4 5]"), "bj");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [nb nc nd nf nk] should be nonnegative integer number "), "bj");
		error(errmsg);
	end

	nb = n(1); nc = n(2); nd = n(3); nf = n(4);
	
	if (size(n,"*") == 4) then
		nk = 1
	else
		nk = n(5);
	end

    //  storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = z(:,2);
    NDATA = size(UDATA,"*");
    function e = G(p,m)
        e = YDATA - _objfunbj(UDATA,p,nd,nc,nf,nb,nk);
    endfunction
    tempSum = nb+nc+nd+nf
    p0 = linspace(0.5,0.9,tempSum)';
    [var,errl] = lsqrsolve(p0,G,size(UDATA,"*"));
    
    err = (norm(errl)^2);
    opt_err = err;
	resid = G(var,[]);
    b = poly([repmat(0,nk,1);var(1:nb)]',"q","coeff");
    c = poly([1; var(nb+1:nb+nc)]',"q","coeff");
    d = poly([1; var(nb+nc+1:nb+nc+nd)]',"q","coeff");
    f = poly([1; var(nb+nd+nc+1:nd+nc+nf+nb)]',"q","coeff");
    t = idpoly(1,coeff(b),coeff(c),coeff(d),coeff(f),Ts)
    
    //  estimating the other parameters
    [temp1,temp2,temp3] = predict(z,t)
    [temp11,temp22,temp33] = pe(z,t)
    
    estData = calModelPara(temp1,temp11,sum(n(1:4)))
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
endfunction

function yhat = _objfunbj(UDATA,x,nd,nc,nf,nb,nk)
    x=x(:)
     q = poly(0,'q')
    tempSum = nb+nc+nd+nf
    //  making polynomials
    b = poly([repmat(0,nk,1);x(1:nb)]',"q","coeff");
    c = poly([1; x(nb+1:nb+nc)]',"q","coeff");
    d = poly([1; x(nb+nc+1:nb+nc+nd)]',"q","coeff");
    f = poly([1; x(nb+nd+nc+1:nd+nc+nf+nb)]',"q","coeff");
    bd = coeff(b*d); cf = coeff(c*f); fc_d = coeff(f*(c-d));
    if size(bd,"*") == 1 then
        bd = repmat(0,nb+nd+1,1)
    end
    maxDelay = max([length(bd) length(cf) length(fc_d)])
    yhat = [YDATA(1:maxDelay)]
    for k=maxDelay+1:size(UDATA,"*")
        bdadd = 0
        for i = 1:size(bd,"*")
            bdadd = bdadd + bd(i)*UDATA(k-i+1)
        end
        
        fc_dadd = 0
        for i = 1:size(fc_d,"*")
            fc_dadd = fc_dadd + fc_d(i)*YDATA(k-i+1)
        end
        cfadd = 0
        for i = 2:size(cf,"*")
            cfadd = cfadd + cf(i)*yhat(k-i+1)
        end
        yhat = [yhat; [ bdadd + fc_dadd - cfadd ]];
    end
endfunction
