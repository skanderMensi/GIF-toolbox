function [] = mainFullFit()
clear,close('all'),clc

addpath /tools

cName = 'GIF1';
disp(cName);

generateData(cName);
disp('reference dataset: OK');
mainFitSubV(cName);
disp('subthreshold dynamics: OK');
mainFitVThreshold(cName);
disp('threshold dynamics: OK');
mainTestMd(cName);
disp('Md*: OK');
mainVoltageTraces(cName);
 
fileName = sprintf('GIFTraining_%s',cName);
load(fileName)
tempError = sprintf('Error Param = %.2f per cent',GIF.eParam);
tempMd = sprintf('Md* = %.2f',GIF.MdStarGIF);
tempRMSE = sprintf('RMSE(V) = %.2f mV',mean(GIF.rmseTest));
tempCPUTime = sprintf('CPU time = %.2f s',GIF.CPUTime+GIF.CPUTimeMd+GIF.threshold.CPUTime);

disp('GIF Fitted:')
disp(tempError)
disp(tempRMSE)
disp(tempMd)
disp(tempCPUTime)

end