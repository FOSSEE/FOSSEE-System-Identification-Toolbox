function frdplot(varargin)

// Plot the frequency domain response
// 
// Calling Sequence
// frdplot(plantData)
// Parameters
// frdData : frd type module
// 
// Description
// plantData must be frd type. Function takes the frequency and response data, and plot the bode plot.  
// Examples
// frdData = (1:100)';
// respData = rand(100,1) + %i * rand(100,1);
// plantData = frd(frdData,respData);
// frdplot(plantData)
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar


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
