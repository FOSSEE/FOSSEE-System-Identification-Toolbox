
function sys = armaX(varargin)
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

	na = n(1); nb = n(2); nc = n(3); //nd = n(4);nf = n(5);
	
	if (size(n,"*") == 3) then
		nk = 1
	else
		nk = n(4);
	end

    // storing U(k) , y(k) and n data in UDATA,YDATA and NDATA respectively 
    YDATA = z(:,1);
    UDATA = z(:,2);

    sys = estpoly(z,[na,nb,nc,0,0,nk])
    
endfunction
