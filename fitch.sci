function varargout = fitch(varargin)
    [lhs , rhs] = argn()
	if ( rhs <> 1 ) then
		errmsg = msprintf(gettext("%s: Wrong number of input arguments"), "fitch");
		error(errmsg)
	elseif typeof(varargin(1)) <> "idpoly" then
        error(msprintf(gettext("%s:Input model must be ""idpoly"" type.\n"),"fitch"))
    end
    model = varargin(1)
    MSE = model.Report.Fit.MSE
    FPE = model.Report.Fit.FPE
    FitPer = model.Report.Fit.FitPer
    AIC = model.Report.Fit.AIC
    AICc = model.Report.Fit.AICc
    nAIC = model.Report.Fit.nAIC
    BIC = model.Report.Fit.BIC
    t = tlist(['fitch','MSE','FPE','FitPer','AIC','AICc','nAIC','BIC'],MSE,FPE,FitPer,AIC,AICc,nAIC,BIC)
    varargout(1) = t
endfunction

function %fitch_p(mytlist)
    f = fieldnames(mytlist)
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
        objectData = mytlist(f(ii))
        if ceil(objectData)-objectData then
            mprintf("%.4f",objectData)
        else
            mprintf("%d",objectData)
        end
        mprintf("\n")
    end
endfunction
