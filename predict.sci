//[y,x0] = predict(data,idpoly,k)
//References
//Digital Control(11.1.2) by Kanna M.Moudgalya 
//System Identification Theory for User Second Edition (3.2) by Lennart Ljung
// Code Aurthor - Ashutosh Kumar Bhargava

function varargout = predict(varargin)
    [lhs,rhs] = argn(0)
//------------------------------------------------------------------------------
// checking the number of inputs
    if rhs < 2 || rhs > 3 then
        error(msprintf(gettext("%s:Wrong number of input arguments.\n"),"predict"))
    end
//------------------------------------------------------------------------------
    data = varargin(1)
    model = varargin(2)
    if rhs == 3 then
        kStep = varargin(3)
    elseif rhs == 2 then
        kStep = 1
    end
    
//------------------------------------------------------------------------------

// k step analysis
    if typeof(kStep) <> 'constant' || isnan(kStep)  then
        error(msprintf(gettext("%s:Prediction horizon(k) must be a non-negative integer number or inf.\n"),"predict"))
    end
// if given k step is infinity or []
    if isinf(kStep) || ~size(kStep,'*') then
        kStep = 1
    end
// checking the dimensions
    if size(kStep,'*') <> 1 || (ceil(kStep)-kStep) then
        error(msprintf(gettext("%s:Prediction horizon(k) must be a non-negative integer number or inf.\n"),"predict"))
    end
//------------------------------------------------------------------------------
// checking the plant model
    if typeof(model) ~= 'idpoly' then
        error(msprintf(gettext("%s:Plant model must be ""idpoly"" type.\n"),"predict"))
    end
    modelSampleTime = model.Ts
    modelTimeUnit = model.TimeUnit
//------------------------------------------------------------------------------
//checking the data type
    if typeof(data) <> 'iddata' && typeof(data) <> 'constant' then
        error(msprintf(gettext("%s:Sample data must be ""iddata"" type or ""n x 2"" matrix type.\n"),"predict"))
    end
// checking the plant data
    if typeof(data) == 'iddata' then
        if ~size(data.OutputData,'*') || ~size(data.InputData,'*') then
            error(msprintf(gettext("%s:Number of sample data in input and output must be equal.\n"),"predict"))
        end
        plantSampleTime = data.Ts
        plantTimeUnit = data.TimeUnit
        data = [data.OutputData data.InputData]
        //disp('iddata')
    elseif typeof(data) == 'constant' then
        if size(data,'c') ~= 2 then
            error(msprintf(gettext("%s:Number of sample data in input and output must be equal.\n"),"predict"))
        end
        plantSampleTime = model.Ts
        plantTimeUnit = model.TimeUnit
    end
//------------------------------------------------------------------------------
// comparing the sampling time
    if modelSampleTime-plantSampleTime <> 0 then
        error(msprintf(gettext("%s:The sample time of the model and plant data must be equal.\n"),"predict"))
    end
// Comparing the time units
    if ~strcmp(modelTimeUnit,plantTimeUnit)  then
    else
        error(msprintf(gettext("%s:Time unit of the model and plant data must be equal.\n"),"predict"))
    end
//------------------------------------------------------------------------------
// ckecking the k step size. if it greater than number of sample size then the 
// k step will become 1 
    if kStep >= size(data,'r') then
        kStep = 1
    end
//------------------------------------------------------------------------------
    //storing the plant data
    //           B(z)               C(z)
    // y(n) = ---------- u(n) + ---------- e(n)
    //         A(z)*F(z)         A(z)*D(z)
    aPoly = poly(model.a,'q','coeff');
    bPoly = poly(model.b,'q','coeff');
    cPoly = poly(model.c,'q','coeff');
    dPoly = poly(model.d,'q','coeff');
    fPoly = poly(model.f,'q','coeff');
    Gq = bPoly/(aPoly*fPoly)
    Hq = cPoly/(aPoly*dPoly)
    
    if kStep == 1 then
        Wkq = Hq^-1 
    elseif kStep > 1 then
        adCoeff = coeff(aPoly*dPoly);adCoeff = adCoeff(length(adCoeff):-1:1);
        adPoly = poly(adCoeff,'q','coeff')
        cCoeff = model.c;cCoeff = cCoeff(length(cCoeff):-1:1);
        cPoly = poly(cCoeff,'q','coeff')
        hBar = clean((ldiv(cPoly,adPoly,kStep))',0.00001)
        hBarPoly = poly(hBar,'q','coeff')
        Wkq = hBarPoly*Hq^-1
    end
    WkqGq = Wkq * Gq
    tempWkqGq = coeff(WkqGq.den)
    if tempWkqGq(1) <> 1 then
        WkqGq.num = WkqGq.num/tempWkqGq(1)
        WkqGq.den = WkqGq.den/tempWkqGq(1)
    end
    Wkq1 = 1-Wkq
    tempWkq1 = coeff(Wkq1.den)
    if tempWkq1(1) == 1 then
        Wkq1.num = Wkq1.num/tempWkq1(1)
        Wkq1.den = Wkq1.den/tempWkq1(1)
    end
   //pause
//------------------------------------------------------------------------------
    // storing the plant data
    uCoeff = coeff(WkqGq.num*Wkq1.den)
    yCoeff = coeff(WkqGq.den*Wkq1.num)
    yCapCoeff = coeff(WkqGq.den*Wkq1.den)
    //pause
    lengthuCoeff = length(uCoeff)
    lengthyCoeff = length(yCoeff)
    lengthyCapCoeff = length(yCapCoeff)
//------------------------------------------------------------------------------
    // keeping initial conditions equal to zero
    uData = zeros(lengthuCoeff,1)
    yData = zeros(lengthyCoeff,1)
    yCapData = zeros(lengthyCapCoeff-1,1)
    uData = [uData;data(:,2)]
    yData = [yData;data(:,1)]
    sampleData = size(data,'r')
    //pause
    // reversing the coefficients
    if ~size(uCoeff,'*') then
        uCoeff = 0
    else
        uCoeff = uCoeff(lengthuCoeff:-1:1)
    end
    if ~size(yCoeff) then
        yCoeff = 0
    else
        yCoeff = yCoeff(lengthyCoeff:-1:1)
    end
    if ~size(yCapCoeff,'*') then
        yCapCoeff = 0
    else
        yCapCoeff = -yCapCoeff(lengthyCapCoeff:-1:2)
    end
    //pause
    for ii = 1:sampleData+1
         //pause
         if ~size(uData(ii:ii+lengthuCoeff-1),'*') then
             tempu = 0
         else
            tempu = uCoeff*uData(ii:ii+lengthuCoeff-1);
         end
         if ~size(yData(ii:ii+lengthyCoeff-1),'*')
             tempy = 0
         else
            tempy = yCoeff*yData(ii:ii+lengthyCoeff-1);
         end
         if ~size(yCapData(ii:ii+lengthyCapCoeff-2),'*') then
            tempyCap = 0
         else
            tempyCap = yCapCoeff*yCapData(ii:ii+lengthyCapCoeff-2);
         end
         yCapData = [yCapData;tempu+tempy+tempyCap];
    end
   // pause
    extraSample = abs(size(yCapData,'r')-sampleData)
    yCapData = yCapData(extraSample+1:$)
    timeData = ((modelSampleTime:modelSampleTime:(size(yCapData,'r')*modelSampleTime))');
    if lhs == 1 then
        clf()
        plot(timeData,yCapData)
        axisData = gca()
        tempTimeUnit = 'Time('+modelTimeUnit+')'
        xtitle('Predicted Response',tempTimeUnit,'y')
        xgrid
        varargout(1) = 0
    elseif lhs == 2 then
        varargout(1) = yCapData
        varargout(2) = 0
    elseif lhs == 3 then
        varargout(1) = yCapData
        varargout(2) = timeData
        varargout(3) = 0
    end
endfunction
