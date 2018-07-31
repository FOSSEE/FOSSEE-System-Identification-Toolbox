function varargout = frd(varargin)


// Stores frequency and response data
// 
// Calling Sequence
// plantData = frd(respData,frdData,Ts)
// Parameters
// frdData : nx1 matrix of non-decreasing frequency points
// respData : nx1 matrix of the frequency response
// Ts : non-negative real number
// plantData : frd type module
// Description
// It is a frd type module that stores the frequency and response data with sampling time Ts. The time unite is in second and frequency unit rad/sec. 
// Examples
// frdData = 0:1024 
// frdData=frdData';
// respData = rand(1024,1)
// Ts = 0.1
//  plantData = frd(frdData,respData,Ts)
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar


    [lhs,rhs] = argn(0)
    if rhs < 2 || rhs > 4 then
        errmsg = msprintf(gettext("%s: Wrong numbers of input arguments."), "frd");
        error(errmsg)
    end
    frequency = varargin(1)
    freqUnit = 'rad/TimeUnit'
    if ~iscolumn(frequency) then
        errmsg = msprintf(gettext("%s: frequency must be a finite column vector."), "frd");
        error(errmsg)
    end
    respData = varargin(2)
//     pause
    if size(frequency,'r') <> size(respData,'r') then
        errmsg = msprintf(gettext("%s: input output matrix dimension must be equal."), "frd");
        error(errmsg)
    end
    if rhs == 2 then
        Ts = 0
    elseif rhs >2 then
        Ts = varargin(3)
    end
    if Ts < 0 || size(Ts,'*') <> 1 || typeof(Ts) <> 'constant' then
        errmsg = msprintf(gettext("%s: Sampling time must be a scalar non negative real number."), "frd");
        error(errmsg)
    end
    //  saving the spectrum value
    if rhs == 4 then
        spect = varargin(4)
    else
        spect = []
    end
    // / matching its dimensions
    if ~size(spect) then
    elseif size(frequency,'r') <> size(spect,'r') then
        errmsg = msprintf(gettext("%s: Numbers of power spectra must be equal to the numbers of frequency."), "frd");
        error(errmsg)
    end
    TUnit = 'seconds'
    t = tlist(['frd','Frequency','FrequencyUnit','ResponseData','Ts','TimeUnit','Spect'],frequency,freqUnit,respData,Ts,TUnit,spect)
    varargout(1) = t
endfunction

// overloading
function %frd_p(varargin)
    myTlist = varargin(1)
    f = fieldnames(myTlist)
    freqData = myTlist.Frequency
    tempRespData= myTlist.ResponseData
    for jj = 1:size(tempRespData,'c')
        respData = tempRespData(:,jj)
        mprintf("\t -------------------------")
        mprintf("\n")
        mprintf("\t Frequency \t  Response")
        mprintf("\n")
        mprintf("\t -------------------------")
        mprintf("\n")
        for ii = 1:length(myTlist.Frequency)
            temp = ''
            if real(respData(ii))>=0 then
                temp = temp + ' '
            end
            temp = temp + string(real(respData(ii)))
    //         temp = string(real(respData(ii)))
            if imag(respData(ii)) > 0 then
                temp = temp +"+"
            end
            if ~imag(respData(ii)) then
            else
                temp = temp + string(imag(respData(ii))) +"i"
            end
    //         temp = temp + string(imag(respData(ii))) + " i"
            mprintf("\n\t %f \t %s",freqData(ii),temp)// real(respData(ii)),imag(respData(ii)))
        end
        mprintf("\n\n")
    end
    if ~myTlist.Ts then
        mprintf("\n  Continuous Domain frequency response.")
    else
        mprintf("\n  Sampling Time = "+string(myTlist.Ts)+" "+myTlist.TimeUnit)
        mprintf("\n  Discrete Domain frequency response.")
    end
endfunction
