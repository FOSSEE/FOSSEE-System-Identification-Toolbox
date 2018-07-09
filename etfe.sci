function varargout = etfe(varargin)
    [lhs,rhs] = argn()
    data = varargin(1)
    if rhs == 1 then
        n = 128
    elseif rhs == 2 then
        n = varargin(2)
    end
    y = data.OutputData;
    u = data.InputData
    N = size(y,'r')
    v = linspace(1,N,n)
    y1 = y(v);u1 = u(v)
    //y1 = y((1:(N-1)/(n-1):N));u1 = u((1:(N-1)/(n-1):N))
    //y1($) = y(N);u1($) = u(N)
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
    output = frd(frq,resp,1)
    varargout(1)= output
endfunction
