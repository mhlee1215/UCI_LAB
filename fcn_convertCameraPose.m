function [] = fcn_convertCameraPose( dataRoot, klgPath)
%FCN_CONVERTGROUNDNORMAL Summary of this function goes here
%   Detailed explanation goes here

% g = gpuDevice(1);
% reset(g);

disp(sprintf('dataRoot:%s', dataRoot));
disp(sprintf('inPath:%s', klgPath));
% disp(sprintf('outPath:%s', outPath));

% addpath(genpath('/home/mhlee/lab_codes'));
addpath(genpath('/home/mhlee/lab_codes/libs/DML/Quaternions'));
% progressbar2(0);
% maxStep = 5;
% curStep = 1;

cvtInfoFilePath =        sprintf('%s/%s_cvt.ply.info', dataRoot, klgPath);
cvtCameraPoseFilePath =  sprintf('%s/%s_cvt.freiburg', dataRoot, klgPath);
dataCameraPoseFilePath = sprintf('%s/%s.freiburg', dataRoot, klgPath);

R = eye(3);
zOffset = 0;
if exist(cvtInfoFilePath, 'file') == 2    
    infoFileID = fopen(cvtInfoFilePath, 'r');
    C = textscan(infoFileID,'%f');
    pData = C{1};
    R = q_getRotationMatrix(pData(1:4));
    zOffset = pData(5);
end

if exist(dataCameraPoseFilePath, 'file') == 2    
    
    cvtCameraPoseFileId = fopen(cvtCameraPoseFilePath, 'w');
    
    dataCameraPoseFileID = fopen(dataCameraPoseFilePath);
    C = textscan(dataCameraPoseFileID,'%f');
    pData = C{1};
%     fclose(dataCameraPoseFileID);
% %     celldisp(C)
% 
%     poseSet{data_idx} = [];
%     cameraPosMat{data_idx} = [];
%     cameraPosMat2{data_idx} = [];
    frameSize = size(pData, 1)/8;
    for fi=0:frameSize-1
        if mod(fi, 100) == 0
            disp(sprintf('converting.. %d/%d', fi, frameSize-1));
        end
        timeStamp = pData(fi*8+1);
        t = pData(fi*8+2:fi*8+4);
        quanternion = pData(fi*8+5:fi*8+8);
        quant_reorder = [quanternion(4) ;quanternion(1:3)];
        r = q_getRotationMatrix(quant_reorder);
        
        t2 = R * t - [0 0 zOffset]';
        r2 = R * r;
        q2 = q_getFromRotationMatrix(r2);
        
        fprintf(cvtCameraPoseFileId, '%f %f %f %f %f %f %f %f\n', ...
            timeStamp, t2(1), t2(2), t2(3), q2(2), q2(3), q2(4), q2(1));
%         pData(fi*8+5:fi*8+8)
    end
end


fclose(infoFileID);
fclose(dataCameraPoseFileID);
fclose(cvtCameraPoseFileId);
% dataCameraPoseFilePath = sprintf('%s/%s.klg.freiburg', dataRoot, inPath);
% if exist(dataCameraPoseFilePath, 'file') == 2    
%     dataCameraPoseFileID = fopen(dataCameraPoseFilePath);
%     C = textscan(dataCameraPoseFileID,'%f');
%     pData = C{1};
%     fclose(dataCameraPoseFileID);
% %     celldisp(C)
% 
%     poseSet{data_idx} = [];
%     cameraPosMat{data_idx} = [];
%     cameraPosMat2{data_idx} = [];
%     frameSize = size(pData, 1)/8;
%     for fi=0:frameSize-1
%         t = pData(fi*8+2:fi*8+4);
% %         pData(fi*8+5:fi*8+8)
% 
%         quanternion = pData(fi*8+5:fi*8+8);
%         quant_reorder = [quanternion(4) ;quanternion(1:3)];
%         r = q_getRotationMatrix(quant_reorder);
%         cameraLoc = -r'*t;
%         cameraLoc2 = t;
% 
%         poseSet{data_idx}(end+1).t = t;
%         poseSet{data_idx}(end).r = r;
%         poseSet{data_idx}(end).cameraLoc = cameraLoc;
%         cameraPosMat{data_idx} = [cameraPosMat{data_idx} cameraLoc];
%         cameraPosMat2{data_idx} = [cameraPosMat2{data_idx} cameraLoc2];
%     end
% 
% end
%     
    
disp('Finished.');

end

