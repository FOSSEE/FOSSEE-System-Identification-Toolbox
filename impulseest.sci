function varargout = impulseest(varargin)
    [lhs,rhs] = argn(0)
    //checking the number of inputs
    if rhs > 2  then
        error(msprintf(gettext("%s: Unexpected number of input arguments "),"impulseest"))
    end
    modelData = varargin(1)
    if typeof(modelData) <> "idpoly" then
        error(msprintf(gettext("%s: Plant model must be ""idpoly"" type. "),"impulseest"))
    end
    //adding noise
    if rhs == 2 then
        noiseFlag = varargin(2)
        if typeof(noiseFlag) <> 'boolean' then
            error(msprintf(gettext("%s: Last input data must be ""boolean"".type "),"impulseest"))
        end
    else
        noiseFlag = %F
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
    if size(find(gsort(abs(roots(sys.den)))>1),'*') then
        tempRoot = clean(roots(sys.den))
        [sorted index] = gsort(abs(real(tempRoot)))
        tempRoot = tempRoot(index)
        n = 26/log10(abs(real(tempRoot(1))))
         
    else
        finalValue = horner(sys,1)*1
        n = 10
        tempValue = sum(ldiv(sys.num,sys.den,n))
        if finalValue >= 0 then
            while finalValue*0.99 >= tempValue || n >=1000
                n = n+5
                tempValue = sum(ldiv(sys.num,sys.den,n))
            end
        elseif finalValue < 0 then
            while finalValue*0.99 <tempValue || n >=1000
                n = n+5
                tempValue = sum(ldiv(sys.num,sys.den,n))
            end
            if n > 100 then
                n = 100
            end  
        end
    end
    uData = [1 zeros(1,n)]
    yData = flts(uData,sys)
    timeData = (0:(n))*modelData.Ts
//    pause
    if noiseFlag then
        varargout(1) = yData'
    else
        stem(timeData,yData)
        h = gcf()
        
        h.figure_name= "Plant Impulse Response"
        varargout(1) = []
    end
endfunction
