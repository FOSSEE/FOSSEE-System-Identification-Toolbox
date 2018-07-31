
function varargout = compare(varargin)
    
// Make comparision between plant output and estimated output
// 
// Calling Sequence
// compare(plantData,sys)
// [yData,tData,fData] = compare(plantData,sys)
// Parameters
// plantData : iddata type or nx2 matrix
// sys : idpoly type polynomial 
// yData : one step ahead estimated output response
// tData : time series data
// fData : real number represent the normalized root mean square (NRMS) of the estimated output response
// Description
// compare function do the estimation of the one step ahead output response of the iddata type polynomial and compare it with the actual output response data with (NRMS). 
// Examples
//  a = [1 0.2];b = [0 0.2 0.3];
//  sys = idpoly(a,b,'Ts',0.1)
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  y = sim(u,sys)+rand(1024,1)
//  plantData = iddata(y,u,0.1)
//  compare(plantData,sys)
// Authors
// Ashutosh Kumar Bhargava 

    [lhs,rhs] = argn(0)
    sampleData = varargin(1)
    sysData = varargin(2)
    yData = []
    uData = []
    Ts = 0
    if lhs > 3 then
        error(msprintf(gettext("%s:Wrong number of output arguments.\n"),"compare"))
    end
    if type(varargin(1)) == 0 then
        error(msprintf(gettext("%s:Null sample data matrix.\n"),"compare"))
    end
    //  if the plant data stores in "iddata" function
    if typeof(sampleData) == 'iddata' then
        yData = sampleData.OutputData
        uData = sampleData.InputData
        Ts = sampleData.Ts
        
    //  if the plant data stores in a matrix of nx2 -->([yData,uData]) 
    elseif typeof(sampleData) == 'constant' then
        if size(sampleData,'c') ~= 2 then
            error(msprintf(gettext("%s:Incorrect number of columns in sample data.\n"),"compare"))
        end
        yData = sampleData(:,1)
        uData = sampleData(:,2)
        Ts = 1
    else
        error(msprintf(gettext("%s:Improper sample datas.\n"),"compare"))
    end
    
    if ~size(yData,'*') then
        error(msprintf(gettext("%s:Sample data must contain plant output datas.\n"),"compare"))
    end
    
    if ~size(uData,'*') then
        error(msprintf(gettext("%s:Sample data must contain plant input datas.\n"),"compare"))
    end
    
    if typeof(sysData) == 'state-space' then
    elseif typeof(sysData) == 'rational' then
        sysData = tf2ss(sysData)
    elseif typeof(sysData) == 'idpoly' then
    else
        error(msprintf(gettext("%s: Wrong type for input argument \n#%d: State-space \n#%d: Transfer function \n#%d: System identification \n expected.\n"),"compare",1,2,3))
    end
    sampleLength = size(yData,'r')
    x0 = []// initial state
    // for state space type systems
    if typeof(sysData) == 'state-space' then
//         sampleLength = size(yData,'r')
        sysData = dscr(sysData,Ts)
        ySys = [0]
        for ii = 2:sampleLength
            tempData = 0
            for jj = 1:ii-1
                tempData = tempData + (sysData.c)*(sysData.a)^(ii-jj-1)*(sysData.b)*uData(jj)
            end
            ySys = [ySys; tempData]
        end
        x0 = zeros(size(sysData.A,'r'),1)
    
    //  for system identification type system(OE,BJ system Models)
    elseif typeof(sysData) == 'idpoly' then
        Ts = sysData.Ts
        zd = [yData uData]
        ySys = compareBJ(sysData,zd)
        // oe model comparision
//         if size(sysData.a,'*') == 1 & size(sysData.c,'*') == 1 & size(sysData.d,'*') == 1 & size(sysData.b,'*') ~= 1 & size(sysData.f,'*') ~= 1 then
//             ySys = compareOE(sysData,zd)
//         else
//             ySys = compareBJ(sysData,zd)
//         end 
    end
    if isrow(yData) then
        yData = yData'
    end
    if isrow(ySys) then
        ySys = ySys'
    end
    fit = fitValue(yData,ySys)
    
    if lhs == 1  then
        varargout(1) = fit// zeros(size(sysData.A,'r'),1)
        dataTime = 0:Ts:(sampleLength-1)*Ts
        plot(dataTime',ySys)
        plot(dataTime',yData,'m')
        // plot(ySys)
        // plot(yData,'m')
        xData = 'Time('
        if typeof(sampleData) == 'constant' then
            xData = xData + 'seconds)'
        elseif typeof(sampleData) == 'iddata' then
            xData = xData + sampleData.TimeUnit+')'
        end
        h = gcf()
        h.figure_name = 'Simulated Compare Plot'
        xgrid()
        xtitle('Comparison of Time Response',xData,'Amplitude')
        plant  = 'plant : ' + string(fit)+ '% fit'
        legend(['sys',plant])
    elseif lhs == 2 then
        varargout(2) = fit
        varargout(1) = x0
    elseif lhs == 3 then
        varargout(3) = fit// zeros(size(sysData.A,'r'),1)
        varargout(2) = x0
        yOutput = iddata(ySys,[],Ts)
        if typeof(varargin(1)) == 'constant' then
            yOutput.TimeUnit = 'seconds'
        elseif typeof(varargin(1)) == 'iddata' then
            yOutout.TimeUnit = sampleData.TimeUnit
        end
        varargout(1) = yOutput
    end
    
endfunction

