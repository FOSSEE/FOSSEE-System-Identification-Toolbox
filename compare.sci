function varargout = compare(varargin)
    [lhs,rhs] = argn(0)
//------------------------------------------------------------------------------
// checking the number of inputs
    if rhs < 2 || rhs > 3 then
        error(msprintf(gettext("%s:Wrong number of input arguments.\n"),"compare"))
    end
//------------------------------------------------------------------------------
    data = varargin(1)
    model = varargin(2)
    if rhs == 3 then
        kStep = varargin(3)
    elseif rhs == 2 then
        kStep = %inf
    end
    
//------------------------------------------------------------------------------

// k step analysis
    if typeof(kStep) <> 'constant' || isnan(kStep)  then
        error(msprintf(gettext("%s:Prediction horizon(k) must be a non-negative integer number or inf.\n"),"compare"))
    end
// if given k step is infinity or []
//    if isinf(kStep) || ~size(kStep,'*') then
//        kStep = 1
//    end
// checking the dimensions
//    if size(kStep,'*') <> 1 || (ceil(kStep)-kStep) then
//        error(msprintf(gettext("%s:Prediction horizon(k) must be a non-negative integer number or inf.\n"),"predict"))
//    end
//------------------------------------------------------------------------------
// checking the plant model
    if typeof(model) ~= 'idpoly' then
        error(msprintf(gettext("%s:Plant model must be ""idpoly"" type.\n"),"compare"))
    end
    modelSampleTime = model.Ts
    modelTimeUnit = model.TimeUnit
//------------------------------------------------------------------------------
//checking the data type
    if typeof(data) <> 'iddata' && typeof(data) <> 'constant' then
        error(msprintf(gettext("%s:Sample data must be ""iddata"" type or ""n x 2"" matrix type.\n"),"compare"))
    end
// checking the plant data
    if typeof(data) == 'iddata' then
        if ~size(data.OutputData,'*') || ~size(data.InputData,'*') then
            error(msprintf(gettext("%s:Number of sample data in input and output must be equal.\n"),"compare"))
        end
        plantSampleTime = data.Ts
        plantTimeUnit = data.TimeUnit
        data = [data.OutputData data.InputData]
        //disp('iddata')
    elseif typeof(data) == 'constant' then
        if size(data,'c') ~= 2 then
            error(msprintf(gettext("%s:Number of sample data in input and output must be equal.\n"),"compare"))
        end
        plantSampleTime = model.Ts
        plantTimeUnit = model.TimeUnit
    end
//------------------------------------------------------------------------------
// comparing the sampling time
    if modelSampleTime-plantSampleTime <> 0 then
        error(msprintf(gettext("%s:The sample time of the model and plant data must be equal.\n"),"compare"))
    end
// Comparing the time units
    if ~strcmp(modelTimeUnit,plantTimeUnit)  then
    else
        error(msprintf(gettext("%s:Time unit of the model and plant data must be equal.\n"),"compare"))
    end
//------------------------------------------------------------------------------
// ckecking the k step size. if it greater than number of sample size then the 
// k step will become 1 
//    if kStep >= size(data,'r') then
//        kStep = 1
//    end
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
    //disp(kStep)
    if kStep == %inf then
        //disp('in inf')
        outData = sim(data(:,2),model)
    else
        if kStep == 1 then
            Wkq = Hq^-1
        elseif kStep > 1 then
            adCoeff = coeff(aPoly*dPoly);adCoeff = adCoeff(length(adCoeff):-1:1);
            adPoly = poly(adCoeff,'q','coeff')
            cCoeff = model.c;cCoeff = cCoeff(length(cCoeff):-1:1);
            cPoly = poly(cCoeff,'q','coeff')
            hBar = clean((ldiv(cPoly,adPoly,kStep))',0.00001)
            hBarPoly = poly(hBar,'q','coeff')
            Wkq = hBarPoly*Hq^-1//*bPoly/(dPoly*fPoly)
        end
    //    pause
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
        z = poly(0,'z');
        WkqGqNum = horner(WkqGq.num,1/z);WkqGqDen = horner(WkqGq.den,1/z);
        uPoly = WkqGqNum/WkqGqDen;
        uPoly = syslin('d',uPoly);uData = flts(data(:,2)',uPoly);
        
        Wkq1Num = horner(Wkq1.num,1/z);Wkq1Den = horner(Wkq1.den,1/z);
        yPoly = Wkq1Num/Wkq1Den;
        yPoly = syslin('d',yPoly);yData = flts(data(:,1)',yPoly);
        outData = (uData+yData)'
    end

    tData = (1:size(data,'r'))'*plantSampleTime
    fitData = fitValue(data(:,1),outData)
    if lhs == 1 then
        varargout(1) = []
        plot(tData,data(:,1),'m')
        plot(tData,outData,'b')
        legend('Plant Data','Model Data : '+string(fitData)+'%')
        xtitle('Comparison of Time Response','Time ('+ plantTimeUnit+')','Amplitude')
        xgrid()
    elseif lhs == 2 then
        varargout(1) = fitData
        varargout(2) = outData
    elseif lhs == 3 then
        varargout(1) = outData
        varargout(2) = tData
        varargout(3) = fitData
    end
    
endfunction
