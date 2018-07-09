function varargout = rarx(varargin)
    [lhs,rhs] = argn(0)
//    
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
    //tempData = [zeros(1,na+nb);tempData]
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
