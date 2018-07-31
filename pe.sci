function varargout = pe(varargin)

// K-steps ahead prediction error
// 
// Calling Sequence
// pe(plantData,sys)
// pe(plantData,sys,k)
// [yData,tData,fData] = pe(plantData,sys)
// [yData,tData,fData] = pe(plantData,sys,k)
// Parameters
// plantData : iddata type or nx2 matrix
// sys : idpoly type polynomial
// k : non-neagtive integer prediction step,default value is 1
// yData : k step ahead prediction error
// tData : time series data
// fData : initial state
// Description
// pe function estimate the differencr between k step ahead predicted and the actual output response. Initial conditions of the dynamic model is zero. 
// Examples
// a = [1 0.2];b = [0 0.2 0.3];
// sys = idpoly(a,b,'Ts',0.1)
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// y = sim(u,sys)+rand(1024,1)
// plantData = iddata(y,u,0.1)
// pe(plantData,sys)
// figure();clf();
// k = 5
// pe(plantData,sys,k)
// Authors
// Ashutosh Kumar Bhargava  

    [lhs,rhs] = argn(0)
    
    
// ------------------------------------------------------------------------------
    data = varargin(1)
    model = varargin(2)
    if rhs == 3 then
        kStep = varargin(3)
    elseif rhs == 2 then
        kStep = 1
    end
// ------------------------------------------------------------------------------
// k step analysis
    if typeof(kStep) <> 'constant' || isnan(kStep)  then
        error(msprintf(gettext("%s:Prediction horizon(k) must be a non-negative integer number or inf.\n"),"pe"))
    end
// if given k step is infinity or []
    if isinf(kStep) || ~size(kStep,'*') then
        kStep = 1
    end
// checking the dimensions
    if size(kStep,'*') <> 1 || (ceil(kStep)-kStep) then
        error(msprintf(gettext("%s:Prediction horizon(k) must be a non-negative integer number or inf.\n"),"pe"))
    end
// ------------------------------------------------------------------------------
// checking the plant model
    if typeof(model) ~= 'idpoly' then
        error(msprintf(gettext("%s:Plant model must be ""idpoly"" type.\n"),"pe"))
    end
    modelSampleTime = model.Ts
    modelTimeUnit = model.TimeUnit
// ------------------------------------------------------------------------------
// checking the data type
    if typeof(data) <> 'iddata' && typeof(data) <> 'constant' then
        error(msprintf(gettext("%s:Sample data must be ""iddata"" type or ""n x 2"" matrix type.\n"),"pe"))
    end
// checking the plant data
    if typeof(data) == 'iddata' then
        if ~size(data.OutputData,'*') || ~size(data.InputData,'*') then
            error(msprintf(gettext("%s:Number of sample data in input and output must be equal.\n"),"pe"))
        end
        plantSampleTime = data.Ts
        plantTimeUnit = data.TimeUnit
        data = [data.OutputData data.InputData]
        // disp('iddata')
    elseif typeof(data) == 'constant' then
        if size(data,'c') ~= 2 then
            error(msprintf(gettext("%s:Number of sample data in input and output must be equal.\n"),"pe"))
        end
        plantSampleTime = model.Ts
        plantTimeUnit = model.TimeUnit
    end
// ------------------------------------------------------------------------------
// comparing the sampling time
    if modelSampleTime-plantSampleTime <> 0 then
        error(msprintf(gettext("%s:The sample time of the model and plant data must be equal.\n"),"pe"))
    end
// Comparing the time units
    if ~strcmp(modelTimeUnit,plantTimeUnit)  then
    else
        error(msprintf(gettext("%s:Time unit of the model and plant data must be equal.\n"),"pe"))
    end
// ------------------------------------------------------------------------------
    sampleLength = size(data,'r')
    // one step ahead prediction
    [y1 ,x0] = predict(data,model)
    errorData1 = data(:,1)-y1
    if kStep == 1 then
        eCapData = errorData1
    elseif kStep > 1 then
        aPoly = poly(model.a,'q','coeff');
        dPoly = poly(model.d,'q','coeff');
        adCoeff = coeff(aPoly*dPoly);adCoeff = adCoeff(length(adCoeff):-1:1);
        adPoly = poly(adCoeff,'q','coeff')
        cCoeff = model.c;cCoeff = cCoeff(length(cCoeff):-1:1);
        cPoly = poly(cCoeff,'q','coeff')
        hBar = clean((ldiv(cPoly,adPoly,kStep))',0.00001)
        hBar = hBar(length(hBar):-1:1)
        hBarLength = length(hBar)
        errorData1 = [zeros(hBarLength,1);errorData1]
        eCapData = []
        
        for ii = 1:sampleLength+1
            eCapData = [eCapData; hBar*errorData1(ii:ii+hBarLength-1)]
        end
    end
    
    timeData = (modelSampleTime:modelSampleTime:(sampleLength)*modelSampleTime)'
    pseudoData = size(eCapData,'r')
    eCapData = eCapData(abs(pseudoData-sampleLength)+1:$)
    // pause
    if lhs == 1 then
        clf()
        plot(timeData,eCapData)
        axisData = gca()
        
        tempTimeUnit = 'Time('+modelTimeUnit+')'
        xtitle('Predicted Response',tempTimeUnit,'y')
        xgrid
        // pause
        varargout(1) = 0
    
    elseif lhs == 2 then
        varargout(1) = eCapData
        varargout(2) = 0
    elseif lhs == 3 then
        varargout(1) = eCapData
        varargout(2) = timeData
        varargout(3) = 0
    end
endfunction
