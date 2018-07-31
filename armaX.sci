
function sys = armaX(varargin)
// Parameters Estimation of ARMAX model using Input Output time-domain data
// 
// Calling Sequence
// sys = armaX(ioData,[na nb nc nk])
// 
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// na : non-negative integer number specified as order of the polynomial A(z^-1)
// nb : non-negative integer number specified as order of the polynomial B(z^-1)+1
// nc : non-negative integer number specified as order of the polynomial C(z^-1)
// nk : non-negative integer number specified as input output delay, Default value is 1
// sys : idpoly type polynomial have estimated coefficients of A(z^-1),B(z^-1) and C(z^-1) polynomials
// 
// Description
// Fit ARMAX model on given input output data 
// The mathematical equation of the ARMAX model 
// <latex>   
// begin{eqnarray}
// A(q)y(n) = B(q)u(n) + C(q)e(n)
// end{eqnarray}
// </latex>
// It is SISO type model. It minimizes the sum of the squares of nonlinear functions using Levenberg-Marquardt algorithm.
// 
// sys ,an idpoly type class, have different fields that contains estimated coefficients, sampling time, time unit and other estimated data in Report object.
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.5];b = [0 2 3];
//  model = iddata(a,b,'Ts',0.1)
//  y = sim(u,model) + rand(length(u),1)
//  plantData = iddata(y,u,0.1)
//  sys = armaX(plantData,[2,2,1])
// 
// Examples
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  a = [1 0.5];b = [0 2 3];
//  model = iddata(a,b,'Ts',0.1)
//  y = sim(u,model) + rand(length(u),1)
//  plantData = [y,u]
//  sys = armaX(plantData,[2,2,1])
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar  

	[lhs , rhs] = argn();	
	if ( rhs < 2 ) then
			errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2"), "armaX", rhs);
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
		errmsg = msprintf(gettext("%s: input and output data matrix should be of size (number of data)*2"), "armaX");
		error(errmsg);
	end

	if (~isreal(z)) then
		errmsg = msprintf(gettext("%s: input and output data matrix should be a real matrix"), "armaX");
		error(errmsg);
	end

	n = varargin(2)
	if (size(n,"*")<3| size(n,"*")>4) then
		errmsg = msprintf(gettext("%s: The order and delay matrix [na nb nc nk] should be of size [3 or 4]"), "armaX");
		error(errmsg);
	end

	if (size(find(n<0),"*") | size(find(((n-floor(n))<%eps)== %f))) then
		errmsg = msprintf(gettext("%s: values of order and delay matrix [na nb nc nk] should be nonnegative integer number "), "armaX");
		error(errmsg);
	end

	na = n(1); nb = n(2); nc = n(3); // nd = n(4);nf = n(5);
	
	if (size(n,"*") == 3) then
		nk = 1
	else
		nk = n(4);
	end

    //  storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = z(:,2);

    sys = estpoly(z,[na,nb,nc,0,0,nk])
    
endfunction
