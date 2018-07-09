function []=Get(varargin)
    //disp('Get')
    tempData = varargin($)
    f = fieldnames(tempData)
//    disp(size(f,'*'))
    if ~size(f,'*') then
        if size(tempData,'*') == 0 then
            printString = 'empty'
        elseif isnan(tempData) then
            printString = 'Nan'
        elseif isinf(tempData) then
            printString = 'inf'
        else
            printString = typeof(tempData)
        end
        mprintf('The type of ""Input"" is : %s \n',printString)
        error(msprintf(gettext("%s: There is no objects in the ""Input"".\n"),"Get"))
    else
        maxLength = []
        for ii = 1:size(f,'*')
            maxLength = [maxLength length(f(ii))]
        end
        maxLength = max(maxLength)
        for ii = 1:size(f,'*')
            blanckSpace = ' '
            for jj = 1:maxLength-length(f(ii))
                blanckSpace = blanckSpace + ' '
            end
            mprintf('\t%s%s : ',blanckSpace,f(ii))
            objectData = tempData(f(ii))
            if typeof(objectData) == 'ce' then
                //mprintf('cell')
                sizeData = size(objectData)
                prodData = prod(sizeData)
                //pause
                if prodData == 1 then
                    objectData = (cell2mat(objectData))
                    objectData = getString(objectData)
                    tempString = '{' + objectData + '}'
                    mprintf('%s',tempString)
                else
                    objectData = getString(sizeData,[])
                    objectData = strsubst(objectData,']',' cell ]')
                    mprintf('%s',objectData)
                end
            elseif typeof(objectData) == 'string' then
                if length(objectData) == 0 then
                    mprintf("{''''}")
                else
                    mprintf('%s ',objectData)
                end
                
            elseif typeof(objectData) == 'constant' then
                //disp('in constant')
                //
                if ~size(objectData,'*') then
                    mprintf('%s%s','[',']')
                    //pause
                    //mprintf('%d',objectData)
                elseif size(objectData,'*') == 1 then
                    if round(objectData)-objectData == 0 then
                        mprintf('%d',objectData)
                    else
                        mprintf('%.2f',objectData)
                    end
                else
                    sizeData = size(objectData)
                    objectData = getString(sizeData,[])
                    objectData = strsubst(objectData,']',' double ]')
                    mprintf('%s',objectData)
                end
            elseif typeof(objectData) == 'polynomial' then
                if size(objectData,'*') == 1 then
                    polyData = pol2str(objectData)
                    mprintf('%s',polyData)
                else
                    objectData = size(objectData)
                    objectData = getString(objectData,[])
                    objectData = strsubst(objectData,']',' polynomial ]')
                    mprintf('%s',objectData)
                end
            elseif typeof(objectData) == 'hypermat' then
                    objectData = size(objectData)
                    objectData = getString(objectData,[])
                    objectData = strsubst(objectData,']',' hypermat  ]')
                    mprintf('%s',objectData)
            end
            mprintf('\n')
        end
    end
endfunction


// converting from matrix to string matrix
function Output = getString(varargin)
    [lhs,rhs] = argn(0)
    if rhs == 1 then
        inBetween = ' '
    elseif rhs == 2 then
        inBetween = 'x'
    end
    stringData = '['
    matData = varargin(1)
    for ii = 1 : length(matData)
        stringData = stringData + string(matData(ii)) //+ inBetween
        if ii == length(matData) then
        else
            stringData = stringData + inBetween
        end
    end
    stringData = stringData + ']'
    Output = stringData
endfunction
