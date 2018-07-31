function sys = idpoly(varargin)

// Stores the identification coefficients
// 
// Calling Sequence
// sys = idpoly(aCoeff,bCoeff,cCoeff,dCoeff,fCoeff,Ts)
// 
// Parameters
// aCoeff : 1xn coefficient matrix of the polynomial A(z^-1)
// bCoeff : 1xn coefficient matrix of the polynomial B(z^-1)
// cCoeff : 1xn coefficient matrix of the polynomial C(z^-1)
// dCoeff : 1xn coefficient matrix of the polynomial D(z^-1)
// fCoeff : 1xn coefficient matrix of the polynomial F(z^-1)
// 
// Description
// It is a idpoly type polynomials that stores the identification coefficients A,B,C,D, and F with sampling time Ts. The time unite is in second. It stores other parameters
//  in Report object.
// Examples
// aCoeff = [1 rand(1,2)]
// bCoeff = [0 rand(1,3)]
// cCoeff = [1 rand(1,2)]
// dCoeff = [1 rand(1,2)]
// fCoeff = [1 rand(1,2)]
// Ts = 0.1
// sys1 = idpoly(aCoeff,bCoeff,cCoeff,dCoeff,fCoeff,Ts)
// sys2 = idpoly(aCoeff,bCoeff,cCoeff,'Ts',Ts)
// sys1 = idpoly(1,bCoeff,cCoeff,1,fCoeff,Ts)
// 
// Authors
// Ashutosh Kumar Bhargava  

    [lhs,rhs] = argn(0)
    tempCell = cell(6,1)
    tempTs = 0
    for ii = 1:rhs
        if typeof(varargin(ii)) == 'string' & varargin(ii) == 'Ts' then
            tempCell{6,1} = varargin(ii+1)
            break
        elseif ~size(varargin(ii),'*') then
            tempCell{ii,1} = 1
        elseif size(varargin(ii),'*') then
            tempCell{ii,1} = varargin(ii)
        end
    end
    // disp(tempCell)
    for ii = 1:6 
        if ~size(cell2mat(tempCell(ii,1)),'*') & ii ~= 6 then
            tempCell{ii,1} = 1
        elseif ~size(cell2mat(tempCell(ii,1)),'*') & ii == 6 then
            tempCell{ii,1} = -1
        end
    end
    // storing the data in A,B,C,D,F matrix
    //           B(z)             D(z)
    //  y(n)=---------- u(n)+ ----------- e(n)
    //        A(z)F(z)          A(z)D(z)
    A = cell2mat(tempCell(1,1));  B = cell2mat(tempCell(2,1));
    C = cell2mat(tempCell(3,1));  D = cell2mat(tempCell(4,1));
    F = cell2mat(tempCell(5,1)); Ts = cell2mat(tempCell(6,1));
    if A(1,1) ~= 1 then
        error(msprintf(gettext("%s: The first coefficient of ""A(z)"" polynomial must be 1.\n"),"idpoly"))
    elseif C(1,1) ~= 1 then
        error(msprintf(gettext("%s: The first coefficient of ""C(z)"" polynomial must be 1.\n"),"idpoly"))
    elseif D(1,1) ~= 1 then
        error(msprintf(gettext("%s: The first coefficient of ""D(z)"" polynomial must be 1.\n"),"idpoly"))
    elseif F(1,1) ~= 1 then
        error(msprintf(gettext("%s: The first coefficient of ""F(z)"" polynomial must be 1.\n"),"idpoly"))
    end
    report = struct('MSE',0,'FPE',0,'FitPer',0,'AIC',0,'AICc',0,'nAIC',0,'BIC',0)
    errors = [0 0 0]
    report = struct('Fit',report,'Uncertainty',errors)
    t = tlist(['idpoly','a','b','c','d','f','Variable','TimeUnit','Ts','Report'],A,B,C,D,F,'z^-1','seconds',Ts,report)
    // t = tlist(['idpoly','a','b','c','d','f','Variable','TimeUnit','Ts'],A,B,C,D,F,'z^-1','seconds',Ts)
    
    sys = t
endfunction


function %idpoly_p(mytlist)
    f = fieldnames(mytlist)
    // A polynomial
    if mytlist(f(1)) == 1 && size(mytlist(f(1)),'*') == 1 then
    else
        mprintf('\n  A(z) =')
        temp = poly2str(mytlist(f(1)))
        mprintf('%s\n',temp)
    end
    
    // B polynomial
    if mytlist(f(2)) == 1 then
    else
        mprintf('\n  B(z) =')
        temp = poly2str(mytlist(f(2)))
        mprintf('%s\n',temp)
    end
    
    // C polynomial
    if mytlist(f(3)) == 1 && size(mytlist(f(3)),'*') == 1 then
    else
        mprintf('\n  C(z) =')
        temp = poly2str(mytlist(f(3)))
        mprintf('%s\n',temp)
    end
    // D polynomial
    if mytlist(f(4)) == 1 && size(mytlist(f(4)),'*') == 1 then
    elseif size(mytlist(f(4)),'*') > 1 then
        mprintf('\n  D(z) =')
        temp = poly2str(mytlist(f(4)))
        mprintf('%s\n',temp)
    end
    
    // F polynomial
    if mytlist(f(5)) == 1 && size(mytlist(f(5)),'*') == 1 then
    else
        mprintf('\n  F(z) =')
        temp = poly2str(mytlist(f(5)))
        mprintf('%s\n',temp)
    end
    
    mprintf('\n  Sampling Time = ')
    
    if mytlist.Ts == -1 then
        mprintf('undefined')
    else
        if (ceil(mytlist.Ts)-mytlist.Ts) == 0 then
            mprintf('%d %s',mytlist.Ts,mytlist.TimeUnit)
        else
            mprintf('%0.4f %s',mytlist.Ts,mytlist.TimeUnit)
        end
    end
    // disp(mytlist.Ts)
    mprintf('\n')
    if mytlist.Report.Fit.MSE then
        temp = ['MSE','FPE','FitPer','AIC','AICc','nAIC','BIC']
        spaces = ' '
        for ii = 1:size(temp,'c')
            digiLength = length(string(round(mytlist.Report.Fit(temp(ii)))))
            digiLength = digiLength + 5-length(temp(ii))
            blank = ''
            for jj = 1:digiLength+1
                blank = blank + " "
            end
            spaces = spaces+blank+temp(ii)+' '
        end
        mprintf('\n')
        mprintf(spaces)
        mprintf("\n  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f",mytlist.Report.Fit.MSE,mytlist.Report.Fit.FPE,mytlist.Report.Fit.FitPer,mytlist.Report.Fit.AIC,mytlist.Report.Fit.AICc,mytlist.Report.Fit.nAIC,mytlist.Report.Fit.BIC)
    end
    
        
endfunction


function strout = poly2str(h)
    temp = poly(h,'x','coeff')
    temp = pol2str(temp)
    temp = strsubst(temp,'-',' - ')
    temp = strsubst(temp,'x^',' z^-')
    temp = strsubst(temp,'x',' z^-1')
    temp = strsubst(temp,'*','')
    temp = strsubst(temp,'+',' + ')
    [ind which]=strindex(temp,'-')
    // disp(ind)
//     disp(which)
    if ind(1,1) ~= 2 then
        temp = ' ' + temp
    elseif ind(1,1) == 2 then
        temp = part(temp,4:length(temp))
        temp = ' -' + temp
    end
    strout = temp
endfunction
