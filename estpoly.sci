

function sys = estpoly(varargin)

// Parameters Estimation of generalized discrete time model using Input Output time-domain data
// 
// Calling Sequence
// sys = bj(ioData,[na nb nc nd nf nk])
// 
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// na : non-negative integer number specified as order of the polynomial A(z^-1)
// nb : non-negative integer number specified as order of the polynomial B(z^-1)+1
// nc : non-negative integer number specified as order of the polynomial C(z^-1)
// nd : non-negative integer number specified as order of the polynomial D(z^-1)
// nf : non-negative integer number specified as order of the polynomial F(z^-1)
// nk : non-negative integer number specified as input output delay, Default value is 1
// sys : idpoly type polynomial have estimated coefficients of A(z^-1), B(z^-1),C(z^-1),D(z^-1) and F(z^-1) polynomials
// 
// Description
// Fit generalized discrete time model on given input output data 
// The mathematical equation of the generalized discrete time model 
// <latex>   
// begin{eqnarray}
// A(q)y(n) = \frac {B(q)}{D(q)}u(n) + \frac {C(q)}{D(q)}e(t)
// end{eqnarray}
// </latex>
// It is SISO type model. It minimizes the sum of the squares of nonlinear functions using Levenberg-Marquardt algorithm.
// 
// sys ,an idpoly type class, have different fields that contains estimated coefficients, sampling time, time unit and other estimated data in Report object.
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.2];b = [0 0.2 0.3];
//  model = idpoly(a,b,'Ts',0.1)
//  y = sim(u,model) + rand(length(u),1)
//  ioData = iddata(y,u,0.1)
//  sys = estpoly(ioData,[2,2,2,2,2,1])
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.2];b = [0 0.2 0.3];
//  model = idpoly(a,b,'Ts',0.1)
//  y = sim(u,model) + rand(length(u),1)
//  ioData = [y,u]
//  sys = estpoly(ioData,[2,2,2,2,2,1])
// 
// Authors
// Ashutosh Kumar Bhargava

	[lhs , rhs] = argn();	
	if ( rhs < 2 ) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "estpoly", rhs);
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
		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "estpoly");
		error(errmsg);
	end

	if (~isreal(z)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "estpoly");
		error(errmsg);
	end

	n = varargin(2)
	if (size(n,"*")<5| size(n,"*")>6) then
		errmsg = msprintf(gettext("%s: The order and delay matrix [na nb nc nd nf nk] should be of size [5 6]"), "estpoly");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [nb nc nd nf nk] should be nonnegative integer number "), "estpoly");
		error(errmsg);
	end

	na = n(1); nb = n(2); nc = n(3); nd = n(4);nf = n(5);
	
	if (size(n,"*") == 5) then
		nk = 1
	else
		nk = n(6);
	end

    //  storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = z(:,2);
    NDATA = size(UDATA,"*");
    function e = G(p,m)
        e = YDATA - _oestpolyfun(UDATA,p,na,nb,nc,nd,nf,nk);// _oestpolyfun(UDATA,p,nd,nc,nf,nb,nk);
    endfunction
    tempSum = na+nb+nc+nd+nf
    p0 = linspace(0.0001,0.001,tempSum)';
    [var,errl] = lsqrsolve(p0,G,size(UDATA,"*"));
    // disp(errl)    
    err = (norm(errl)^2);
    // disp(err)
    opt_err = err;
	resid = G(var,[]);
    x = var
//     b = poly([repmat(0,nk,1);var(1:nb)]',"q","coeff");
//     c = poly([1; var(nb+1:nb+nc)]',"q","coeff");
//     d = poly([1; var(nb+nc+1:nb+nc+nd)]',"q","coeff");
//     f = poly([1; var(nb+nd+nc+1:nd+nc+nf+nb)]',"q","coeff");
    a = poly([1; x(1:na)]',"q","coeff");
    b = poly([repmat(0,nk,1);x(na+1:na+nb)]',"q","coeff");
    c = poly([1; x(na+nb+1:na+nb+nc)]',"q","coeff");
    d = poly([1; x(na+nb+nc+1:na+nb+nc+nd)]',"q","coeff");
    f = poly([1; x(na+nb+nd+nc+1:na+nd+nc+nf+nb)]',"q","coeff");
    t = idpoly(coeff(a),coeff(b),coeff(c),coeff(d),coeff(f),Ts)
    
        // t = sys;// idpoly(1,coeff(b),coeff(c),coeff(d),coeff(f),Ts)
    
    //  estimating the other parameters
    [temp1,temp2,temp3] = predict(z,t)
    [temp11,temp22,temp33] = pe(z,t)
    // pause
    estData = calModelPara(temp1,temp11,na+nb+nc+nd+nf)
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
    
    // sys.TimeUnit = unit
endfunction

function yhat = _oestpolyfun(UDATA,x,na,nb,nc,nd,nf,nk)// (UDATA,x,nd,nc,nf,nb,nk)
    x=x(:)
    q = poly(0,'q')
    tempSum = na+nb+nc+nd+nf
    //  making polynomials
    a = poly([1; x(1:na)]',"q","coeff");
    b = poly([repmat(0,nk,1);x(na+1:na+nb)]',"q","coeff");
    c = poly([1; x(na+nb+1:na+nb+nc)]',"q","coeff");
    d = poly([1; x(na+nb+nc+1:na+nb+nc+nd)]',"q","coeff");
    f = poly([1; x(na+nb+nd+nc+1:na+nd+nc+nf+nb)]',"q","coeff");
    bd = coeff(b*d); cf = coeff(c*f); fc_d = coeff(f*(c-a*d));
    if size(bd,"*") == 1 then
        bd = repmat(0,nb+nd+1,1)
    end
    // pause
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
        // pause
        yhat = [yhat; [ bdadd + fc_dadd - cfadd ]];
    end
endfunction
