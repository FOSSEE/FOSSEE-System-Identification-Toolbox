function varargout = misdata(varargin)
    [lhs,rhs] = argn(0)
//------------------------------------------------------------------------------
// checking the number of inputs
    if rhs <> 1 then
        error(msprintf(gettext("%s:Wrong number of input arguments.\n"),"misdata"))
    end
//------------------------------------------------------------------------------
    ioData = varargin(1)
    if typeof(ioData) <> "iddata" then
        error(msprintf(gettext("%s: Plant input data must be ""iddata"" type. "),"misdata"))
    end
    inputMat = ioData.InputData;inputMat = linearINTRP(inputMat,abs(ioData.Ts));ioData.InputData = inputMat;
    outputMat = ioData.OutputData;outputMat = linearINTRP(outputMat,abs(ioData.Ts));ioData.OutputData = outputMat;
    varargout(1) = ioData
endfunction

function varargout = linearINTRP(matData,Ts)
    // looking for overall nan values
    nanData = isnan(matData);nanIndex = find(nanData == %T)
    if ~size(nanIndex,'*') then
        varargout(1) = matData
    else
        tempMat = []
        matSize = size(matData,'r')
        // looking for nan in each column 
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
