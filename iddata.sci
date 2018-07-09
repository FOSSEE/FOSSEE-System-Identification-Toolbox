function varargout = iddata(varargin)
    [lhs,rhs] = argn(0)
    if rhs == 2 | rhs == 1 then
        Ts = 1
    elseif rhs == 3 then 
        Ts = varargin(3)
    else
        error(msprintf(gettext("%s:Incorrect number of input arguments.\n"),"iddata"))
    end
    if Ts <= 0 || typeof(Ts) <> 'constant' || ~size(Ts,'*') then
        error(msprintf(gettext("%s:Inconsist sampling time ""Ts"", ""non negative real number"" expected.\n"),"iddata"))
    end
    if rhs == 1 then
        OutputData = varargin(1)
        InputData = []
    elseif rhs == 2 | rhs == 3 then
        OutputData = varargin(1)
        InputData = varargin(2)
        if isrow(InputData) then
            InputData = InputData'
        end
    end
    if size(OutputData,'*') & size(InputData,'*') then
        if size(OutputData,'r') ~= size(InputData,'r') then
            error(msprintf(gettext("%s:The numbers of the plant out datas must be equal to the numbers of the plant input datas.\n"),"iddata"))
        end
    end
    t = tlist(['iddata','OutputData','InputData','Ts','TimeUnit'],OutputData,InputData,Ts,'seconds')
    varargout(1) = t
endfunction

function %iddata_p(mytlist)
    f = fieldnames(mytlist)
    if ~size(mytlist(f(1)),'*') & ~size(mytlist(f(2)),'*') then
        mprintf('  Empty sample data.\n')
    else
        outputSize = size(mytlist(f(1)))
        inputSize = size(mytlist(f(2)))
        if prod(outputSize) then
            sampleSize = max(outputSize)
        elseif prod(inputSize) then
            sampleSize = max(inputSize)
        end
        mprintf('  Time domain sample data having %d samples.',sampleSize)
        if round(mytlist(f(3)))-mytlist(f(3)) == 0 then
            mprintf('\n  Sampling Time = %d',mytlist(f(3)))
        else
            mprintf('\n  Sampling Time = %f',mytlist(f(3)))
        end
        mprintf(' %s',mytlist(f(4)))
        mprintf('\n')
        if prod(outputSize) then
            mprintf('\n  Output channel \n')
            for ii = 1:min(outputSize) 
                yString = 'y' + string(ii)
                mprintf('  %s\n',yString)
            end 
        end
        if prod(inputSize) then
            mprintf('\n  Input channel \n')
            for ii = 1:min(inputSize) 
                uString = 'u' + string(ii)
                mprintf('  %s\n',uString)
            end 
        end
    end
    mprintf('\n')
endfunction
