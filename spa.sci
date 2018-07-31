function varargout = spa(varargin)
// Use spectral analysis to estimate frequency response  
// 
// Calling Sequence
// frdData = spa(plantData,winSize,Freq)
// Parameters
// plantData : iddata type
// winSize : non-neative integer number
// Freq : frequency points to evaluate the response
// frdData : frd type object
// Description
// spa function does the estimation of the frequency response of iddata of a plant using the spectral analysis. Hanning window is used in spectral 
// analysis. Default size of the window is minimum of the number of sample divided by 10 or 30. The default frquency point is %pi x (1:128)/(sampling time x 128).
// Examples
// a = [1 0.2];b = [0 0.2 0.3];
// sys = idpoly(a,b,'Ts',0.1)
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// y = sim(u,sys)+rand(1024,1)
// plantData = iddata(y,u,0.1)
// frdData = spa(plantData)
// Authors
// Bhushan Manjarekar, Ashutosh Kumar Bhargava

    [lhs , rhs] = argn();
    if ( rhs < 1 || rhs > 3 ) then
        errmsg = msprintf(gettext("%s: Wrong number of input arguments" ),"spa");
        error(errmsg)
    elseif typeof(varargin(1)) <> "iddata" then
        error(msprintf(gettext("%s:Plant data must be ""iddata"" type.\n"),"spa"))
    end
    
    plantData = varargin(1)
    windowSize = %F
    inputFreq = %F
// ------------------------------------------------------------------------------
// arranging the plant data    
    inputData = plantData.InputData; 
    if ~size(inputData,"*") then
        error(msprintf(gettext("%s:Input data must be non-zero vector. \n"),"spa"))
    end
    outputData = plantData.OutputData
    if ~size(outputData,"*") then
        error(msprintf(gettext("%s:Output data must be non-zero vector. \n"),"spa"))
    end
// ------------------------------------------------------------------------------
    N = size(inputData,'r')
    nout = size(outputData,'c')
    nin = size(inputData,'c')
    
    if ~windowSize then 
        windowSize = min(N/10, 30)
    end
    
    if inputFreq then
    else 
        inputFreq = (1:128)/128 * %pi/plantData.Ts
    end
    
    M = windowSize
    Ryu = xcov(outputData,inputData,M,'biased')
    Ruu = xcov(inputData,inputData,M,'biased')
    Ryy = xcov(outputData,outputData,M,'biased')
    G = [];spect = [];
    for ii = 1:nout
        phiY = spaCalculation(inputFreq,Ryy,M)
        temp = phiY
        for jj = 1:nin
            phiYU =  spaCalculation(inputFreq,Ryu,M)// sapply(freq, cov2spec, Ryu[i, j, ], M)
            phiU = spaCalculation(inputFreq,Ruu,M)
            G =  [G phiYU./phiU]
            // pause
            temp = temp - phiYU .* conj(phiYU)./phiU
        end
        spect = [spect temp]
    end
//    disp(size(spect))
//    disp(spect)
    frdData = frd(G',inputFreq',plantData.Ts,spect')
    varargout(1) = frdData
endfunction

function varargout = spaCalculation(varargin)
    inputFreq = varargin(1)
    xcovData = varargin(2)
    Mnumber = varargin(3)
    temp = []
    win_l=window('hn',2*Mnumber+1)
    for ii = 1:size(inputFreq,'c')
        seq1 = exp(-(%i) * (-Mnumber:Mnumber) * inputFreq(ii))
        data = (seq1.*win_l)
        data2 = sum(data.*xcovData')
        temp = [temp data2]
    end
    varargout(1) = temp
endfunction
