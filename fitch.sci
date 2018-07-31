function varargout = fitch(varargin)

// Characteristics of estimated model
// 
// Calling Sequence
// estData = fitch(sys)
// Parameters
// sys : idpoly type polynomial
// estData : stuct type variable have following objects
// MSE : Mean Square Error, show the wellness of the estimated model on plant data
// FPE : Final Prediction Error
// FitPer : Normalized root mean squared error (NRMSE), show the precentage fit of the estimated model on plant data
// AIC : Raw Akaike Information Citeria (AIC), show the model quality
// AICc : Small sample-size corrected AIC
// nAIC : Normalized AIC
// BIC : Bayesian Information Criteria (BIC)
// Description
// fitch function represent other calculated model characteristics to show the wellness of the model on different scales.  
// Examples
//  a = [1 0.2];b = [0 0.2 0.3];
//  sys1 = idpoly(a,b,'Ts',0.1)
//  u = idinput(1024,'PRBS',[0 1/20],[-1 1])
//  y = sim(u,sys)+rand(1024,1)
//  plantData = iddata(y,u,0.1)
//  sys = arx(plantData,[2,2,1])
//  estData = fitch(sys)
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar  

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
