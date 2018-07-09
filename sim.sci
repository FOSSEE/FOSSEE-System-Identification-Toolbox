function varargout = sim(varargin)
    [lhs,rhs] = argn(0)
    //checking the number of inputs
    if rhs < 2 || rhs > 3  then
        error(msprintf(gettext("%s: Unexpected number of input arguments "),"sim"))
    end
    modelData = varargin(2)
    inputData = varargin(1)
    //checking the first input
    if typeof(inputData) <> "constant" && typeof(inputData) <> "iddata"  then
        error(msprintf(gettext("%s: Plant input data must be ""iddata"" type or ""double vector"". "),"sim"))
    end
    if typeof(inputData) == "iddata" then
        inputData = inputData.InputData
    end
    if ~iscolumn(inputData) then
        error(msprintf(gettext("%s: Plant input data must be ""double column vector"". "),"sim"))
    end
    //checking the plant model type
    if typeof(modelData) <> "idpoly" then
        error(msprintf(gettext("%s: Plant model must be ""idpoly"" type. "),"sim"))
    end
    //adding noise
    if rhs == 3 then
        noiseFlag = varargin(3)
        if typeof(noiseFlag) <> 'boolean' then
            error(msprintf(gettext("%s: Last input data must be ""boolean"".type "),"sim"))
        end
    else
        noiseFlag = %F
    end
    //adding noise of zero mean and 1 as standard deviation
    if noiseFlag then
        numberOfSamples = size(inputData,'r')
        R = = grand(numberOfSamples,1,"nor",0,1)
        inputData = inputData+R
    end
    z = poly(0,'z')
    aPoly = poly(modelData.a(length(modelData.a):-1:1),'z','coeff')
    bPoly = poly(modelData.b,'z','coeff')
    fPoly = poly(modelData.f(length(modelData.f):-1:1),'z','coeff')
    afPoly = aPoly*fPoly
    
    bCoeff = modelData.b
    extra = 1
    if ~bCoeff(1,1) then
        afCoeff = coeff(afPoly)
        bLength = length(bCoeff);afLength = length(afCoeff)
        if bLength == afLength then
            extra = 1
        else
            extra = z^-(bLength-afLength)
        end
    end
    
    bPoly = poly(modelData.b(length(modelData.b):-1:1),'z','coeff')
    sys = syslin('d',bPoly,afPoly)*extra
    outputData = (flts(inputData',sys))'
    
    if typeof(varargin(1)) == "iddata" then
        outputData = iddata(outputData,[],varargin(1).Ts)
        outputData.TimeUnit = varargin(1).TimeUnit
    end
    varargout(1) = outputData
endfunction
