function sys = dataSlice(data,Start,End,Freq)
    [lhs,rhs] = argn()
    // storing the model data 
    modelData = data
    // storing the statrting point
    try
        startData = Start
    catch
        startData = 1
    end
    //storing the end point
    try 
        endData = End
    catch
        endData = LastIndex(data)
    end
    //Storing the frequency
    try
        freqData = Freq
    catch
        freqData = 1
    end
    // error message generate
    if startData > endData  then
        error(msprintf(gettext("%s:Start index can not greater than End index.\n"),"dataSlice"))
    end
    if size(startData,'*') ~= 1 then
        error(msprintf(gettext("%s:Start index must be non negative scalar integer number.\n"),"dataSlice"))
    end
    if size(endData,'*') ~= 1 then
        error(msprintf(gettext("%s:End index must be non negative scalar integer number.\n"),"dataSlice"))
    end
    if ~freqData || size(freqData,'*') ~= 1 then
        error(msprintf(gettext("%s:Frequency must be non negative scalar number.\n"),"dataSlice"))
    end
    //--------------------------------------------------------------------------
    if typeof(modelData) == 'constant' then
        Ts = 1
    elseif typeof(modelData) == 'iddata' then
        Ts = modelData.Ts
    end
    //--------------------------------------------------------------------------
    if freqData> Ts || modulo(Ts,freqData) then
        warning(msprintf(gettext("%s:inconsistent frequency.\n"),"dataSlice"))
        freqData = Ts
    end
    if typeof(modelData) == 'constant' then
        temp = modelData(startData:Ts/freqData:endData,:)
    elseif typeof(modelData) == 'iddata' then
        tempY = modelData.OutputData;tempU = modelData.InputData
        tempY = tempY(startData:Ts/freqData:endData,:);tempU = tempU(startData:Ts/freqData:endData,:)
        temp = iddata(tempY,tempU,Ts/freqData)
        temp.TimeUnit = modelData.TimeUnit
    end
    sys = temp
endfunction

function varargout = LastIndex(modelData)
    //finding the sample size
    if typeof(modelData) == "constant" then
        varargout(1) = length(modelData(:,1))
        
    elseif typeof(modelData) == "iddata"  then
        temp = modelData.OutputData
        varargout(1) = length(temp(:,1))
    end    
endfunction
