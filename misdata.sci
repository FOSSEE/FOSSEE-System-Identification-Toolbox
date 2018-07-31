function varargout = misdata(varargin)

// Recover Missing Data by Interpolation
// 
// Calling Sequence
// data = misdata(plantData)
// Parameters
// plantData : iddata type object with missing data as Nan
// data : iddata type object with interpolated data
// 
// Description
// misdata function recovers the experimental missing plant time series data by linear interpolation. 
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// u(100:105) = %nan;u(300:100:1000) = %nan
// y(420:422) = %nan;y(555:100:1024) = %nan
// plantData = iddata(y,u,0.1)
// data = misdata(plantData)
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar  


    [lhs,rhs] = argn(0)
// ------------------------------------------------------------------------------
//  checking the number of inputs
    if rhs <> 1 then
        error(msprintf(gettext("%s:Wrong number of input arguments.\n"),"misdata"))
    end
// ------------------------------------------------------------------------------
    ioData = varargin(1)
    if typeof(ioData) <> "iddata" then
        error(msprintf(gettext("%s: Plant input data must be ""iddata"" type. "),"misdata"))
    end
    inputMat = ioData.InputData;inputMat = linearINTRP(inputMat,abs(ioData.Ts));ioData.InputData = inputMat;
    outputMat = ioData.OutputData;outputMat = linearINTRP(outputMat,abs(ioData.Ts));ioData.OutputData = outputMat;
    varargout(1) = ioData
endfunction

function varargout = linearINTRP(matData,Ts)
    //  looking for overall nan values
    nanData = isnan(matData);nanIndex = find(nanData == %T)
    if ~size(nanIndex,'*') then
        varargout(1) = matData
    else
        tempMat = []
        matSize = size(matData,'r')
        //  looking for nan in each column 
        for ii = 1:size(matData,'c')
            nanData = isnan(matData(:,ii));nanIndex = find(nanData == %T);
            if ~size(nanData,'*') then
                tempMat = [tempMat matData(,ii)]
            else
                timeData = (linspace(1*Ts,matSize*Ts , matSize))';
                nanMat = isnan(matData(:,ii));
                data = matData(:,ii)
                data(nanMat) = interp1(timeData(~nanMat), data(~nanMat), timeData(nanMat));
                tempMat = [tempMat data]
            end
        end
        varargout(1) = tempMat
    end
endfunction
