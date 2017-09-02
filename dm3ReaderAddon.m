% Authors: 
%   Phong Hien Nguyen, MRSEC Summer 2017 REU Fellow
%   Peter Sand
% Version:
%   12
% Last Edit: 
%   09 August 2017

% This program uses dm3ReaderModified.m (modified to write tag files) to 
% extract data from *.dm3 files containing EELS spectra images. Zero-loss 
% peak shifting is performed on high-loss spectra using low-loss spectra 
% data (if provided). Background is subtracted from each high-loss spectra 
% using a user-defined range of values for background curve-fitting. 
% Integration of the spectra is performed using a user defined range
% (LLInt+20). A *.csv file containing the element name is created with a 
% table containing distance values in column 1, and intensity values in 
% column 2. It is created in the same directory as this program file.

% Starting Requirements:
% User-input files are *.dm3 type files only.
% dm3ReaderModified.m and EELSData.m must be included in the same 
% directory as this file.

% 1. Run (F5)
% 2. Specify the settings to be used.
% 2. Select the high-loss spectrum image file
% 3. Select the low-loss spectrum image file (only if user specified that 
% a low-loss spectrum is available. Without a low-loss spectrum image, no 
% zero-loss peak shift can be performed.)

% Primary function
% function [output]=dm3ReaderAddon12()

% Close all windows and clear the workspace
close all; clearvars;

currentFolder=pwd;

% Loading settings
settingNum=5;
useSettings=input('Would you like to load previously saved settings? (Yes = 1, No = 0) \n');
if useSettings==1
    [FileName,PathName]=uigetfile('*.txt','Select the settings you wish to load. ');
    fileIDLoad=fopen([PathName FileName],'r');
    for i=1:1:settingNum
        settings(i)=fscanf(fileIDLoad,'%i',1);
    end
    fclose(fileIDLoad);
    
    % specific settings
    fprintf('\n');
    direction=settings(1);
    fprintf('direction = %i \n',direction);
    LLExist=settings(2);
    fprintf('LLExist = %i \n',LLExist);
    estType=settings(3);
    fprintf('estType = %i \n',estType);
    userDefAxis=settings(4);
    fprintf('userDefAxis = %i \n',userDefAxis);
    norm2Material=settings(5);
    fprintf('norm2Film = %i \n',norm2Material);
elseif useSettings==0
    direction=input('Do you want top-to-bottom information or left-to-right information)? \n(top-to-bottom = 1, left-to-right = 0) \n');
    if direction~=1 && direction~=0
        error('Error: Input must be 1 or 0. ');
    else
    end
    settings(1)=direction;
    
    % specific setting
    LLExist=input('Do you have a low-loss spectrum image? (Yes = 1, No = 0) \n');
    if LLExist~=0 && LLExist~=1
        error('Error: Input must be 1 or 0. ');
    else
    end
    settings(2)=LLExist;

    estType=input('What kind of background estimation function do you want to use? (Disable = 0, Power Law = 1, Exponential = 2) \n');
    if estType~=0 && estType~=1 && estType~=2
        error('Error: Input must be 0, 1, or 2. ');
    else
    end
    settings(3)=estType;

    userDefAxis=input('Do you want to manually create the spacial axis? (Yes = 1, No = 0) \n');
    if userDefAxis~=0 && userDefAxis~=1
        error('Error: Input must be 1 or 0. ');
    else
    end
    settings(4)=userDefAxis;

    norm2Material=input('Do you want to normalize the intensity to those within the film? (Yes = 1, No = 0) \n');
    if norm2Material~=0 && norm2Material~=1
        error('Error: Input must be 1 or 0.');
    else
    end
    settings(5)=norm2Material;

    save(settings,currentFolder);
end

% Generates the general input (and output) structure for high-loss and 
% low-loss EELS spectrum images
HLinput=EELSData; HLoutput=EELSData; LLinput=EELSData;

% Calls dm3Reader.m to extract high-loss (HL) data and assigns output to
% HLinput.Data
fprintf('\n');
fprintf('Please select the high-loss EELS image file. \n');
[HLinput.Data,fileID,fileName]=dm3ReaderCall();

% Save file ID for .csv output file creation
txtStart=length(fileName)-2; txtEnd=length(fileName); fileNameTemp=fileName;
fileNameTemp(1,txtStart:txtEnd)='csv';
csvFileName=fileNameTemp;

% Importing useful tags (HL) and tranposition if desired data is not
% top-to-bottom
fileID=[fileID '_tags.txt'];

if direction==1
    % - Pixel row scale (ladder)
    scaleTag='Scale = '; tagNum=7;
    HLinput.scale(1)=tagFind(fileID,scaleTag,tagNum);
    % -Pixel column scale (horizontal)
    tagNum=8;
    HLinput.scale(2)=tagFind(fileID,scaleTag,tagNum);
    % - eV scale
    tagNum=9;
    HLinput.scale(3)=tagFind(fileID,scaleTag,tagNum);
elseif direction==0
    HLinput=cubeTranspose(HLinput);
    % - Pixel row scale (ladder)
    scaleTag='Scale = '; tagNum=7;
    HLinput.scale(2)=tagFind(fileID,scaleTag,tagNum);
    % -Pixel column scale (horizontal)
    tagNum=8;
    HLinput.scale(1)=tagFind(fileID,scaleTag,tagNum);
    % - eV scale
    tagNum=9;
    HLinput.scale(3)=tagFind(fileID,scaleTag,tagNum);
else
end

%   - eV origin
originTag='Origin = '; tagNum=7;
[HLinput.origin]=tagFind(fileID,originTag,tagNum)*-1;
    
% Calls dm3Reader.m to extract low-loss (LL) data and assigns output to
% LLinput.Data and tranposition if desired data is not top-to-bottom
if LLExist==1
    fprintf('\n');
    fprintf('Please select the low-loss EELS image file. \n');
    [LLinput.Data,fileID]=dm3ReaderCall();
    % Transposition
    if direction==0
        LLinput=cubeTranspose(LLinput);
    else
    end
    % Importing useful tags (LL)
    fileID=[fileID '_tags.txt'];
    % - eV origin
    tagNum=7;
    LLinput.origin=tagFind(fileID,originTag,tagNum)*-1;
    % - eV scale
    tagNum=9;
    LLinput.scale(3)=tagFind(fileID,scaleTag,tagNum);
else
end

fprintf('\n');
elementID=input('Element? Or enter 0 to quit. \n','s'); 
isElement=str2double(elementID);
elementNum=0;
while isElement~=0
    csvFileName=[elementID '_EELS.csv']; 
    elementNum=elementNum+1;

    %   - Make a new copy
    HLoutput=HLinput;
    
    % Step 1: Background removal for each pixel
    if estType==1
        % - Define the range to use for curve-fitting
        fprintf('\nLL and UL are used for every pixel. \n');
        fprintf('Following inputs must be multiples of %.3e. \n',HLinput.scale(3));
        LLEst=input('LL for background removal: ');
        ULEst=input('UL for background removal: ');
        % - Fit the curve and subtract the background
        [preShifteVHL]=eVExtract(HLinput); 
        [HLoutput,Est]=backgroundRemove(HLinput,HLoutput,LLEst,ULEst,preShifteVHL,estType);
    else
    end

    % Step 2: Apply zero-loss peak shift
    %   - Define a unified eV vector
    if LLExist==1
        [preShifteVLL]=eVExtract(LLinput);
        [shiftMatrix,minShift,maxShift]=findShiftMatrix(LLinput,preShifteVLL);
        [HLinput.shift]=shiftMatrix; [HLoutput.shift]=shiftMatrix;
        [unifiedeV]=eVUnify(preShifteVHL,minShift,maxShift);
        [output.unifiedeV]=round(unifiedeV,3);
        % - Reconstruct data cube (HLoutput.Data) to conform to unified eV vector
        for x=1:1:HLoutput.dim(1)
            for y=1:1:HLoutput.dim(2)
                [zTemp]=pad(HLoutput,x,y,minShift,maxShift);
                for z=1:1:length(unifiedeV)
                    HLoutput.Data(x,y,z)=zTemp(z);
                end
            end
        end
    else
    end
    output.Data=HLoutput.Data;
    
    % Step 3: Integration
    %   - Define the range to use for integration
    LLInt=input('LL for integration: ');
    ULInt=LLInt+20; %input('input('UL for integration (multiple of %.3f): \n',HLinput.scale(1));
    %   - Create a structure (intSquare) containing integrated values
    %   (percentages)
    if LLExist==0
        unifiedeV=preShifteVHL;
    else
    end

    LLIntPos=findLocation(unifiedeV,LLInt);
    ULIntPos=findLocation(unifiedeV,ULInt);
    for x=1:1:HLoutput.dim(1)
        for y=1:1:HLoutput.dim(2)
            [output.intSquare(x,y),boundedY]=pixInt(HLoutput,x,y,...
                unifiedeV,LLIntPos,ULIntPos);
        end
        output.intVector=sum(output.intSquare);
        output.intVector=output.intVector./max(output.intVector);
    end

    %   - Creation of pixel array axis
    if userDefAxis==1
        if elementNum<2
            imgHeight=input('What is the physical height of the EELS spectrum image? \n');
            if imgHeight<=0
                error('Error: image must have real height ');
            else
            end
        else
        end
        HLoutput.scale(1)=imgHeight/HLoutput.dim(2);
        HLoutput.scale(2)=HLoutput.scale(1); % pixels are squares, therefore length should equal width
        output.columnLength=lengthExtract(HLoutput,1);
        output.rowLength=lengthExtract(HLoutput,2);
    else
        output.columnLength=lengthExtract(HLoutput,1);
        output.rowLength=lengthExtract(HLoutput,2);
    end

    % Step 4: Normalizing Data
    %   - Creation, processing, and extraction of distance vs integration value 
    % data set. Column 1 contains depth. Column 2 contains corresponding 
    % element intensity.
    output.intVector=[output.rowLength' output.intVector'];
    if norm2Material==1
        if elementNum<2
            figure(1);
            plot(output.intVector(:,2)); xlabel('Pixel'); ylabel('Counts, %');
            xticks([1:1:length(output.intVector(:,2))]); 
            grid on;
            matLLPos=input('At what pixel does the material start? \n');
            matULPos=input('At what pixel does the material end? \n');
            close figure 1;
        else
        end
        output.intVectorNorm=normalize(output.intVector,matLLPos,matULPos);
        output.intVectorNorm=[output.rowLength' output.intVectorNorm];
        % storing user-defined film start and user-defined film end for EDX
        % data comparison
        output.materialCoordinates(1)=matLLPos*HLoutput.scale(2);
        output.materialCoordinates(2)=matULPos*HLoutput.scale(2);
    else
    end

    % Step 5: Extraction
    
    % Extraction Preparation Start-----------------------------------------
    
    % OriginFormat specifies whether or not the BRPixels and dataSquare
    % outputs will be written in left-to-right, top-to-bottom order (0) or
    % right-to-left, bottom-to-top order (1). This is convenient for 
    % plotting in Origin, as applying an offset will show in the more 
    % intuitive left-to-right, top-to-bottom order as you look at the 
    % Origin plot from the top to the bottom. 
    OriginFormat=1;
    
    % format per-pixel background-removed data into an exportable form
    BRPixels=pixelPrint(HLoutput);
    if OriginFormat==1
        BRPixels=squareFlip(BRPixels); 
    else
    end
    BRPixels=[unifiedeV' BRPixels];
    output.BRPixels=BRPixels;
    
    % sum background removed data across rows (top-to-bottom) or columns
    % (left-to-right)
    dataSquare=[rowSum(HLoutput)']; 
    if OriginFormat==1
        dataSquare=squareFlip(dataSquare); 
    else
    end
    dataSquare=[unifiedeV' dataSquare];
    output.dataSquare=dataSquare;
    
    % create a folder for all exported *.csv files
    [status,msg,msgID]=mkdir('dm3ReaderAddon_output');
    if status==1
        cd([pwd '\dm3ReaderAddon_output']);
    else
        error('Folder creation was unsuccessful. Please close all files saved within the folder and try again. \n');
    end
    % Extraction Preparation End-------------------------------------------
    
    % Extract per-pixel background-removed data
    % Column 1 contains the unified eV vector. Remaining columns contain
    % pixel-specific spectra from right-to-left, bottom-to-top order
    % (useful when plotting in Origin and applying a constant offset to all
    % columns plotted. Resulting Origin curves, from top-to-bottom are
    % pixel-specific spectra: in pixel order: left-to-right, top-to-bottom
    % For more intuitive left-to-right, top-to-bottom order, change the
    % OriginFormat value from 1 to 0.
    writeBRPixels=1;
    if writeBRPixels==1
        csvwrite(['BRPixels_' csvFileName],output.BRPixels);
    else
    end
    
    % Extract dataSquare (information summed across rows (top-to-bottom) or
    % columns (left-to-right))
    % Column 1 contains the unified eV vector. Remaining columns contain
    % summed spectra in bottom-to-top or right-to-left order of the EELS 
    % image according to the initial user-specificed input of direction.
    % For more intuitive left-to-right or top-to-bottom order, dchange the 
    % OriginFormat value from 1 to 0.
    writeDataSquare=1;
    if writeDataSquare==1
        csvwrite(['DataSquare_' csvFileName],output.dataSquare);
    else
    end

    % Extract element-specific intensity. If direction = 1, the intensity 
    % of the element in the pixel has the same indices as the pixel in the 
    % EELS image.
    writeIntSquare=1;
    if writeIntSquare==1
        csvwrite(['IntSquare_' csvFileName],output.intSquare);
    else
    end
    
    % Extract normalized-integrated data summed across pixel rows 
    % (top-to-bottom information) or pixel columns (left-to-right 
    % information) or corresponding non-normalized-integrated data 
    % Column 1 contains the unified eV vector.
    writeIntVector=1;
    if norm2Material==1 && writeIntVector==1
        csvwrite(['IntVectorNorm_' csvFileName],output.intVectorNorm);
    elseif norm2Material==0 && writeIntVector==1
        csvwrite(['IntVector_' csvFileName],output.intVector);
    else
        error('writeIntVector must equal 1 (write the normalized element-concentration vs depth information) or 2 (write the non-normalized element-concentration vs depth information.');
    end

    % Plotting
    run=0;
    if run==1 % plots a background-removed spectrum
        % for i=1:1:HLoutput.dim(2)
        figure;
        plot(unifiedeV,squeeze(HLoutput.Data(2,13,:))); % plot(unifiedeV,squeeze(HLoutput.Data(2,i,:)));
        title('Background Removal Test'); 
        xlabel('eV'); ylabel('Counts'); grid on;
        % end

        if norm2Material==1 % plots the element concentration vs depth plot
            figure;
            plot(output.intVectorNorm(:,1),output.intVectorNorm(:,2),'-x');
            title('Intensity vs Distance'); 
            xlabel('Distance, micrometers'); ylabel('y'); grid on;
        else
            figure;
            plot(output.intVector(:,1),output.intVector(:,2),'-x');
            title('Intensity vs Distance'); 
            xlabel('Distance, micrometers'); ylabel('y'); grid on;
        end
    else
    end
    
    cd(currentFolder);
    elementID=input('\n\nNext element? Or 0 to quit. \n','s');
    isElement=str2double(elementID);
end
% end

% -------------------------------------------------------------------------
% Functions Defined Below
% -------------------------------------------------------------------------
% Creates a unified energy vector by padding the left and right sides of
% each pixel's eV spectrum with 0s
function [zTemp]=pad(HLoutput,x,y,minShift,maxShift)
% number of zeros to use on left and right side respectively
    numLeftPad=(HLoutput.shift(x,y)-minShift)/HLoutput.scale(3); %!!!
    numLeftPad=int64(numLeftPad);
    numRightPad=(maxShift-HLoutput.shift(x,y))/HLoutput.scale(3); %!!!
    numRightPad=int64(numRightPad);
% vectors with appropriately indexed zero values
    leftPad=zeros(1,numLeftPad);
    rightPad=zeros(1,numRightPad);
% extracts pixel-specific intensity spectrum
    for i=1:1:HLoutput.dim(3)
        zTemp(i)=HLoutput.Data(x,y,i);
    end
% creates padded pixel-specific intensity spectrum
    zTemp=[leftPad zTemp rightPad];
end

% creates a corresponding unified eV vector to paddedDataCubeHL
% (eV axis is shifted in this function)
function [unifiedeV]=eVUnify(preShifteV,minShift,maxShift)
    startVal=preShifteV(1)-maxShift; % first value
    endVal=preShifteV(length(preShifteV))-minShift; %last value
    stepSize=preShifteV(2)-preShifteV(1); % step size
    unifiedeV=[startVal:stepSize:endVal]; % define output
end

% Generates an eV-axis vector (pre-zero loss peak shift)
function [preShifteV]=eVExtract(inputStructure)
    start=inputStructure.origin(1)*inputStructure.scale(3);
    lasteV=start+(2047*inputStructure.scale(3));
    preShifteV=[start:inputStructure.scale(3):lasteV];
end

% Generates a length vector for pixel row or columns
% n = 1 -> columns
% n = 2 -> rows
function [length]=lengthExtract(inputStructure,n)
   start=inputStructure.scale(n)/2;
   lasteV=start+((inputStructure.dim(n)-1)*inputStructure.scale(n));
   length=[start:inputStructure.scale(n):lasteV];
end

% Finds useful data in the tag file generated by modified dm3Reader.m
% fileID and literal are strings
function [value]=tagFind(fileID,literal,tagNum)
fid=fopen(fileID);
tline=fgetl(fid);
counter=0;
    while ischar(tline)
       matches=strfind(tline,literal);
       num=length(matches);
       if num>0
          counter=counter+1;
          if counter==tagNum
            line=tline;
          else
          end
       else
       end
       tline=fgetl(fid);
    end
fclose(fid);
value=textscan(line,['%s %s %f']);
value=value{3};
end

% Calls dm3Reader
function [dataCube,fileID,fileName]=dm3ReaderCall()
    dm3ReaderModified;
    dataCube=ans.SI;
    fileID=ans.fileID;
    fileName=ans.fileName;
end

function [HLoutput,zApprox]=backgroundRemove(HLinput,HLoutput,LLEst,ULEst,...
    preShifteVHL,estType)
    for x=1:1:HLinput.dim(1)
        for y=1:1:HLinput.dim(2)
            LLPos=findLocation(preShifteVHL,LLEst);
            ULPos=findLocation(preShifteVHL,ULEst);
            xTemp2=preShifteVHL(1,LLPos:1:ULPos);
            for i=LLPos:1:ULPos
                if HLinput.Data(x,y,i)>=1
                    yTemp(i-LLPos+1)=HLinput.Data(x,y,i);
                else
                    yTemp(i)=yTemp(i-1);
                end
            end
            if estType==0 % disable
                dataCubeBR=squeeze(HLinput.Data(x,y,:));
            elseif estType==1 % power law
                p=polyfit(log(xTemp2),log(yTemp),1); m=p(1); b=exp(p(2));
                zData=squeeze(HLinput.Data(x,y,:)); zData=zData';
                for z=1:1:length(zData)
                    zApprox(x,y,z)=b.*preShifteVHL(z).^m;
                    dataCubeBR(z)=zData(z)-zApprox(x,y,z);
                end
            elseif estType==2 % decaying exponential
                p=polyfit(xTemp2,log(yTemp),1); m=p(1); b=exp(p(2));
                zData=squeeze(HLinput.Data(x,y,:)); zData=zData';
                for z=1:1:length(zData)
                    zApprox(x,y,z)=b.*exp(m.*preShifteVHL(z));
                    dataCubeBR(z)=zData(z)-zApprox(x,y,z);
                end
            else
            end
            HLoutput.Data(x,y,:)=dataCubeBR(1,:);
            clear dataCubeBR;
        end
    end       
end

% finds the position of a particular input
function [location]=findLocation(vector,value)
    position=0; % initial counter position
    for i=1:1:length(vector)
        position=position+1;
        temp=vector(i);
        temp=int64(temp);
        if temp==value
            location=position;
        else
        end
    end
end

% Sums the filtered data across rows, creating a data square in dataSquare.
% Intensity values are preserved in columns. Row identity is preserved.
function [dataSquare]=rowSum(HLoutput)
    for y=1:1:HLoutput.dim(2)
        for z=1:1:HLoutput.dim(3)
            zNet(z)=sum(HLoutput.Data(:,y,z));
            dataSquare(y,z)=zNet(z);
        end
    end
end

% Creates a shift-matrix, with each element containing the estimated
% appropriate shift of the corresponding pixel
function [shiftMatrix,minShift,maxShift]=findShiftMatrix(LLinput,preShifteVLL)
    for x=1:1:LLinput.dim(1)
        for y=1:1:LLinput.dim(2)
            shiftMatrix(x,y)=shiftFinder(LLinput,preShifteVLL,x,y);
            maxShift=max(max(shiftMatrix));
            minShift=min(min(shiftMatrix));
        end
    end
end

% Finds the appropriate shift using the estimated position of the greatest
% intensity in the low-loss EELS spectrum
function [shift]=shiftFinder(LLinput,preShifteVLL,x,y)
    for z=1:1:LLinput.dim(3)
        zTemp(z)=LLinput.Data(x,y,z);
    end
    [M,I]=max(zTemp);
    i=1;
    zMaxPos(i)=I; 
    for z=I+1:1:LLinput.dim(3)
        if LLinput.Data(x,y,z)==M
            i=i+1;
            zMaxPos(i)=z;
        else
        end
    end
    shiftPos=round(sum(zMaxPos)/i);
    shift=preShifteVLL(shiftPos);
end      

% Integrates a pixel-specific EELS spectrum given the pixel coordinate
function [intVal,boundedY]=pixInt(HLoutput,x,y,unifiedeV,LLIntPos,ULIntPos)
    boundedX=unifiedeV(1,LLIntPos:ULIntPos);
    boundedY=squeeze(HLoutput.Data(x,y,LLIntPos:ULIntPos));
    intVal=trapz(boundedX,boundedY);
end

% Normalizes the integrated intensities within a specific range
function [intVectorNorm]=normalize(intVector,matLLPos,matULPos)
    intVectorNorm=intVector(:,2)-min(intVector(:,2));
    temp=intVectorNorm(matLLPos:matULPos,1);
    intVectorNorm=intVectorNorm/max(temp);
end

% Transposes cubic datasets (with EELSData-defined fields), 
% switching only the first two indices
function [outputCube]=cubeTranspose(inputCube)
    outputCube=EELSData;
    for x=1:1:inputCube.dim(1)
        for y=1:1:inputCube.dim(2)
            outputCube.Data(y,x,:)=inputCube.Data(x,y,:);
        end
    end
end

% Flip data from bottom-to-top to top-to-bottom (useful for Origin)
% Note: flipping the cube left-to-right flips image data top-to-bottom
function [outputCube]=cubeFlip(inputCube)
for i=1:1:inputCube.dim(2)
    outputCube.Data(:,i,:)=inputCube.Data(:,inputCube.dim(2)-i+1,:);
end
end

% Flip data from bottom-to-top to top-to-bottom (useful for Origin)
% Note: flipping the square left-to-right flips image data top-to-bottom
function [outputSquare]=squareFlip(inputSquare)
tempSize=size(inputSquare);
for i=1:1:tempSize(2)
    outputSquare(:,i)=inputSquare(:,tempSize(2)-i+1);
end
end

% Formats background-subtracted 3D data into 2D data
% column 1 = eV-axis, remaining columns contain pixel-specific EELS, 
% from left-to-right, top-to-bottom
function [BRPixels]=pixelPrint(HLoutput)
BRPixels=[];
for x=1:1:HLoutput.dim(1)
    for y=1:1:HLoutput.dim(2)
        BRPixels=[BRPixels squeeze(HLoutput.Data(x,y,:))];
    end
end
end

% Prompts user for the decision to save settings and the name for the 
% settings file.
function save(settings,currentFolder)
saveSettings=input('Would you like to save these settings? (Yes = 1, No = 0) \n');
if saveSettings==1
    directoryName=uigetdir('C:\','Select a folder in which to save your settings ');
    cd(directoryName);
    settingsName=input('What do you want to name your settings file? (No spaces or special characters) \n','s');
    delete [settingsName '.txt'];
    fileIDPrint=fopen([directoryName '\' settingsName '.txt'],'a');
    fprintf(fileIDPrint,'%i ',settings);
    fclose(fileIDPrint);
else
end
% Change back to working folder
cd(currentFolder);
end
