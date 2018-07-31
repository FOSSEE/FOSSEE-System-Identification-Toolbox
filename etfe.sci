function varargout = etfe(varargin)

// Estimate empirical transfer function
// 
// Calling Sequence
// frdData = etfe(plantData,n)
// Parameters
// plantData : iddata type
// n : frequency sample spacing, default value is 128
// frdData : frd type object
// Description
// etfe function takes time domain plant data,iddata type, and estimate the empirical transfer function by taking the ratio of the fourier transforms
// of the output and the input variables
// Examples
//  a = [1 0.2];b = [0 0.2 0.3];
//  sys = idpoly(a,b,'Ts',0.1)
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  y = sim(u,sys)+rand(1024,1)
//  plantData = iddata(y,u,0.1)
//  frdData = etfe(plantData)
// Authors
// Ashutosh Kumar Bhargava,Bhushan Manjarekar  

    [lhs,rhs] = argn()
    if rhs > 2 || rhs < 1 then
        error(msprintf(gettext("%s:Wrong number of input arguments.\n"),"etfe"))
    end
    data = varargin(1)
    if typeof(data) <> "iddata" then
        error(msprintf(gettext("%s:Plant time series data must be ""iddata"" type.\n"),"etfe"))
    end
    if rhs == 1 then
        n = 128
    elseif rhs == 2 then
        n = varargin(2)
    end
    y = data.OutputData;
    u = data.InputData;
    if ~size(u,'r') then
        error(msprintf(gettext("%s:Non zero input data point needed.\n"),"etfe"))
    elseif ~size(y,'r') then
        error(msprintf(gettext("%s:Non zero output data point needed.\n"),"etfe"))
    end
    N = size(y,'r')
    
    y1 = y((1:ceil((N-1)/(n-1)):N));u1 = u((1:ceil((N-1)/(n-1)):N))
    y1($) = y(N);u1($) = u(N)
    data12 = [y1,u1]
    z = [fft(y1),fft(u1)]
    z = z/size(z,'r')
    magData1 = abs(z(:,1));magData2 = abs(z(:,2))
    argData1 = phasemag(z(:,1),'m');argData2 = phasemag(z(:,2),'m')
    magData = magData1./magData2;argData = argData1-argData2
    argData = [cosd(argData) sind(argData)]
    data = [magData.*argData(:,1) magData.*argData(:,2)]
    output = data(:,1)+%i*data(:,2)
    resp = output(1:ceil(length(output)/2))
    frq = (1: ceil(n/2)) * %pi/floor(n/2)
    output = frd(frq',resp,1)
    varargout(1)= output
endfunction
