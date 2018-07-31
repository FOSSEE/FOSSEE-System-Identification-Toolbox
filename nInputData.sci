function varargout = nInputSeries(varargin)
// Count the number of input output data series
// 
// Calling Sequence
// ioData = nInputSeries(plantData)
// 
// Parameters
// plantData : iddata type module
// ioData : structure type variable having input and output object
// 
// Description
// plantData must be iddata type. The function returns the number of input and output data series in struct form. It has two fields name as input and output.  
// 
// Examples
// plantData = iddata(rand(100,3),rand(100,2),0.1)
// ioData = nInputSeries(plantData)
// 
// Examples
// plantData = iddata(rand(100,2),[],0.1)
// ioData = nInputSeries(plantData)
// 
// Examples
// plantData = iddata([],rand(100,3),0.1)
// ioData = nInputSeries(plantData)
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar
// 
    [lhs rhs] = argn(0)
    if rhs <> 1 then
        error(msprintf(gettext("%s: Wrong number of input arguments."),"nInputSeries"))
    end
    iddataData  =  varargin(1)
    if typeof(iddataData) <> 'iddata' then
        error(msprintf(gettext("%s:Wrong type for input argument %d: ""iddata"" expected."),"nInputSeries",1))
    end
    t = struct('input',0,'output',0)
    
    if ~size(iddataData.InputData,'*') then
        t.input = 0
    else
        t.input = size(iddataData.InputData,'c')
    end
    if ~size(iddataData.OutputData,'*') then
        t.output = 0
    else
        t.output = size(iddataData.InputData,'c')
    end
    varargout(1) = t
endfunction
