function varargout = compareBJ(varargin)
    //varargin(1) -> idpoly data about oe
    //varargin(2) -> [y u] matrix of "nx2" dimension
    //disp('compareBj')
    bjData = varargin(1)
    //disp(typeof(bjData))
    plantData = varargin(2)
    //disp(typeof(plantData))
    yData = plantData(:,1)
    uData = plantData(:,2)
    uData = [0;uData]
    //storing the data in A,B,C,D,F matrix
    //          B(z)             C(z)
    // y(n)=---------- u(n)+ ----------- e(n)
    //       A(z)F(z)          A(z)D(z)
    polyA = poly(bjData.a,'x','coeff')
    polyB = poly(bjData.b,'x','coeff')
    polyC = poly(bjData.c,'x','coeff')
    polyD = poly(bjData.d,'x','coeff')
    polyF = poly(bjData.f,'x','coeff')
    adf = polyA*polyD*polyF
    bd = polyB*polyD
    cf = polyC*polyF
    delay = max(size(coeff(adf),'*'),size(coeff(bd),'*'),size(coeff(cf),'*'))
    yHat = [0]
    bdCoeff = coeff(bd)
    adfCoeff = coeff(adf)
    adfCoeff = -adfCoeff(2:length(adfCoeff))
    for ii = 1:length(uData)
        uSum = 0;ySum = 0;
        for jj = 1:length(bdCoeff)
            if ii-jj <= 0 then
                uSum = uSum + 0
            else
                uSum = uSum + uData(ii-jj+1)*bdCoeff(jj)
            end
        end
        for jj = 1:length(adfCoeff)
            if ii-jj <= 0 then
                ySum = ySum + 0
            else
                ySum = ySum + yHat(ii-jj+1)*adfCoeff(jj)
            end
        end
        yHat = [yHat; uSum+ySum]
    end
    tempStart = 1
    if size(yHat,'r')- size(yData,'r') > 0 then
        tempStart = size(yHat,'r')-size(yData,'r')+1
    end
    varargout(1) = yHat(tempStart:length(yHat))
//    plot(yHat(tempStart:length(yHat)),'m')
//    plot(yData)
//    xgrid()
//    pause
endfunction
