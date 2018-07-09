loadmatfile('rData.mat')
z = [fft(data1(:,1)),fft(data1(:,2))]
z = z/size(z,'r')
magData1 = abs(z(:,1));magData2 = abs(z(:,2))
argData1 = phasemag(z(:,1),'m');argData2 = phasemag(z(:,2),'m')
magData = magData1./magData2;argData = argData1-argData2
argData = [cosd(argData) sind(argData)]
data = [magData.*argData(:,1) magData.*argData(:,2)]
output = data(:,1)+%i*data(:,2)
resp = output(1:ceil(length(output)/2))
frq = (1: ceil(256/2)) * %pi/floor(256/2)
output = frd(frq,resp,1)
disp(output)
