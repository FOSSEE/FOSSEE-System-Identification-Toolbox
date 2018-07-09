function iddataplot(varargin)
    [lhs rhs] = argn(0)
    if rhs <> 1 then
        error(msprintf(gettext("%s: Wrong number of input arguments."),"iddataplot"))
    end
    iddataData  =  varargin(1)
    if typeof(iddataData) <> 'iddata' then
        error(msprintf(gettext("%s:Wrong type for input argument %d: ""iddata"" expected."),"iddataplot",1))
    end
//    figure()
//    xtitle('Plant i/o Data','Time('+iddataData.TimeUnit+')','Amplitude')
    uData = iddataData.InputData
    yData = iddataData.OutputData
    timeLength = max(size(uData,'r'),size(yData,'r'))
    timeData = ((1:1:timeLength)*iddataData.Ts)'
    //timeData = timeData(1:length(iddataData)-1)
    if size(uData,'*') && size(yData,'*') then
        firstIndex = 2
    else
        firstIndex = 1
    end
    // ploting y data
    h = gcf()
    if h.figure_name == "Plant Input Output Data"; then
        clf()
    end
    if size(yData,'*') then
        secondIndex = size(yData,'c')
        if secondIndex == 1 then
            outputString = 'y'
        else
            outputString = []
            for ii = 1:secondIndex
                outputString = [outputString 'y'+string(ii)]
            end
        end
        //disp(outputString)
        for ii = 1:secondIndex
            subplot(firstIndex,secondIndex,ii);plot(timeData,yData(:,ii))
            xtitle(outputString(ii))
        end
    end
    if size(uData,'*') then
        secondIndex = size(uData,'c')
     //   disp(secondIndex)
        if secondIndex == 1 then
            outputString = 'u'
        else
            outputString = []
            for ii = 1:secondIndex
                outputString = [outputString 'u'+string(ii)]
            end
        end
        //disp(outputString)
        for ii = 1:secondIndex
            if size(yData,'*') then
                temp = 1
            else
                temp = 0
            end
            subplot(firstIndex,secondIndex,ii+secondIndex*temp);plot2d(timeData,uData(:,ii),2)
            //pause
            xtitle(outputString(ii))
        end
    end
    h = gcf()
    //disp(h)
    h.figure_name= "Plant Input Output Data";
endfunction
