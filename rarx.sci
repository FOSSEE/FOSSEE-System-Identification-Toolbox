function varargout = rarx(varargin)
    
// Parameters Estimation of ARX model by recursive method
// 
// Calling Sequence
// sys = rarx(ioData,[na nb nk],lambda)
// Parameters
// ioData : iddata or [outputData inputData] ,matrix of nx2 dimensions, type plant data
// na : non-negative integer number specified as order of the polynomial A(z^-1)
// nb : non-negative integer number specified as order of the polynomial B(z^-1)+1
// nk : non-negative integer number specified as input output delay, Default value is 1
// lambda : Forgetting factor,Default value is 0.95
// sys : idpoly type polynomial have estimated coefficients of A(z^-1) and B(z^-1) polynomials
// 
// Description
// Fit RARX model on given input output data 
// RARX model is SISO type model. It uses recursive weighted least-squares algorithm to estimate the coefficient of ARX model 
// sys is a struct type variable output contains data about theta and yhat.
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// ioData = iddata(y,u,0.1)
// sys = rarx(ioData,[2,2,1])
// 
// Examples
// u = idinput(1024,'PRBS',[0 1/20],[-1 1])
// a = [1 0.2];b = [0 0.2 0.3];
// model = idpoly(a,b,'Ts',0.1)
// y = sim(u,model) + rand(length(u),1)
// ioData = [y,u]
// sys = rarx(ioData,[2,2,1])
// 
// Authors
// Ashutosh Kumar Bhargava, Bhushan Manjarekar  

    [lhs,rhs] = argn(0)
    plantData = varargin(1)
    orderData = varargin(2)
    na = orderData(1);nb = orderData(2)
    // arranging na ,nb,nk
    if size(orderData,"*") == 2 then
        nk = 1
    elseif size(orderData,'*') == 3 then
        nk = orderData(3)
    end
    // storing the lambda value
    if rhs == 3 then
        lambda = varargin(3)
    else
        lambda = 0.95
    end
    
    nb1 = nb + nk - 1
    n = max(na, nb1)    
    // arranging the plant data
    if typeof(plantData) == 'constant' then
        Ts = 1;unitData = 'second'
    elseif typeof(plantData) == 'iddata' then
        Ts = plantData.Ts;unitData = plantData.TimeUnit
        plantData = [plantData.OutputData plantData.InputData]
    end
    N = size(plantData,'r')
    uIndex = nk:nb1
    yIndex = []
    if na <> 0 then
        yIndex = 1:na
    end
    df = N - na - nb
    Plast = 10^4 * (eye(na+nb,na + nb))
    theta = zeros(N + 1, na + nb)
    yHat = plantData(:,1);yData = plantData(:,1)
    tempData = zeros(N,na+nb)
    for ii = 1:na
        tempData(ii+1:ii+N,ii) = -plantData(:,1)
    end
    // arranging samples of u matrix
    for ii = 1:nb
        tempData(ii+nk:ii+N+nk-1,ii+na) = plantData(:,2)
    end
    // tempData = [zeros(1,na+nb);tempData]
    tempData = tempData(1:N+1,:)
    
    for ii = 1:N
        temp = tempData(ii,:)
        yHat(ii) = temp*theta(ii,:)'
        eps_i = yData(ii)-yHat(ii)
        kappa_i = Plast * temp'/(lambda + temp * Plast * temp')
        theta(ii+1,:) = ((theta(ii,:))' + eps_i * kappa_i)'
        Plast = (eye(na + nb,na + nb) - kappa_i * (temp)) * Plast/lambda
        
    end
    theta = theta(1:N,:)
    yHat = yHat(1:N)
    
    varargout(1) = struct('theta',theta,'yHat',yHat)
endfunction
