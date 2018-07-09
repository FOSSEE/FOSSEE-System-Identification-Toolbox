function frdplot(varargin)
    [lhs rhs] = argn(0)
    if rhs <> 1 then
        error(msprintf(gettext("%s: Wrong number of input arguments."),"frdplot"))
    end
    frdData  =  varargin(1)
    if typeof(frdData) <> 'frd' then
        error(msprintf(gettext("%s:Wrong type for input argument %d: ""frd"" expected."),"frdplot",1))
    end
    bode((frdData.Frequency)',(frdData.ResponseData)')
endfunction
