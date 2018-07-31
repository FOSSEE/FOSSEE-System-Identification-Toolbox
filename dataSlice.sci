
function sys = dataSlice(data,Start,End,Freq)
// Select sample data from iddata 
// 
// Calling Sequence
// h = dataSlice(plantData,Start,End,Ts)
// Parameters
// data : iddata type  
// Start : non-negative integer index
// End : non-negative integer index, always greater than Start index
// Ts : sampling frequency, default value is 1
// Description
// Extracts the samples in between Start and End index of the plant time series data,iddata type. For specified sampling frequency, it resamples the extracted data.    
// 
// Examples
//  a = [1 0.5];b = [0 0.2 0.3];
//  sys = idpoly(a,b,'Ts',0.1)
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  y = sim(u,sys)+rand(1024,1)
//  plantData = iddata(y,u,0.1)
//  h = dataSlice(plantData,1,20,0.1)
// Authors
// Ashutosh Kumar Bhargava,Bhushan Manjarekar  

    [lhs,rhs] = argn()
    //  storing the model data 
    modelData = data
    //  storing the statrting point
    try
        startData = Start
    catch
        startData = 1
    end
    // storing the end point
    try 
        endData = End
    catch
        endData = LastIndex(data)
    end
    // Storing the frequency
    try
        freqData = Freq
    catch
        freqData = 1
    end
    //  error message generate
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
    // --------------------------------------------------------------------------
    if typeof(modelData) == 'constant' then
        Ts = 1
    elseif typeof(modelData) == 'iddata' then
        Ts = modelData.Ts
    end
    // --------------------------------------------------------------------------
    if freqData> Ts || modulo(Ts,freqData) then
        warning(msprintf(gettext("%s:inconsistent frequency.\n"),"dataSlice"))
        freqData = Ts
    end
    if typeof(modelData) == 'constant' then
        temp = modelData(startData:Ts/freqData:endData,:)
    elseif typeof(modelData) == 'iddata' then
        tempY = modelData.OutputData;tempU = modelData.InputData
        if ~size(tempY,'r') then
            tempY = []
        else
            tempY = tempY(startData:Ts/freqData:endData,:);
        end
        if ~size(tempU,'r') then
            tempU = []
        else
            tempU = tempU(startData:Ts/freqData:endData,:)
        end
        temp = iddata(tempY,tempU,Ts/freqData)
        temp.TimeUnit = modelData.TimeUnit
    end
    sys = temp
endfunction

function varargout = LastIndex(modelData)
    // finding the sample size
    if typeof(modelData) == "constant" then
        varargout(1) = length(modelData(:,1))
        
    elseif typeof(modelData) == "iddata"  then
        temp = max(size(modelData.OutputData,'r'),size(modelData.InputData,'r'))
        varargout(1) = temp
    end    
endfunction
