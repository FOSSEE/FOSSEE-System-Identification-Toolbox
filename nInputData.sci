function varargout = nInputSeries(varargin)
        [lhs rhs] = argn(0)
    if rhs <> 1 then
        error(msprintf(gettext("%s: Wrong number of input arguments."),"iddataplot"))
    end
    iddataData  =  varargin(1)
    if typeof(iddataData) <> 'iddata' then
        error(msprintf(gettext("%s:Wrong type for input argument %d: ""iddata"" expected."),"nInputSeries",1))
    end
    if ~size(iddataData.InputData,'*') then
        varargout(1) = 1
    else
        varargout(1) = size(iddataData.InputData,'c')
    end
endfunction
