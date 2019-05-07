% May need to change these paths if the locations of these modules differ
% on your setup.
addpath('/home/common/matlab/shadedErrorBar');
addpath('/home/common/matlab/colormap');
addpath('../CaCLEAN');
addpath('../CaCLEAN/smoothn');
addpath('../CaCLEAN/FMINSEARCHBND/FMINSEARCHBND');

% Change these for a different model permutation (w/ vs. w/o mitochondria,
% 1um enforced spacing model vs. standard RyR spacing)
filename = 'simulatedMicroscopy_outputs/simulatedMicroscopyResults_interp53_SNR100_1umSpacing_NoMito_54x80x7.mat';
plotTitle = 'Low cluster density (N=408), no mitochondria'; % 984  408
outfile = 'TestCaCLEAN_outputs/LowDensity_NoMito';

% "Figure 3" is an example of classification for visual inspection. Pick a
% slice (between 1-22) and an admissible window (tolerance). Originally set
% with 10nm admissible window resolution (i.e. a tolerance of 64 corresponds 
% with an AW of 640nm)
fig3Slice = 15; 
fig3Tol = 64;

% Load up the data for this example
load(filename);

% Set up some preliminary variables
clusterDiameter = 0.2; % 200 nm

numTol = size(RyrTolerances,2);
numSlices = size(MultiIdenoised,1);
Recall = NaN(numSlices, numTol);
Prec = NaN(numSlices, numTol);

FP_detected = NaN(numSlices, numTol);
FN_actual = NaN(numSlices, numTol);
TP_detected = NaN(numSlices, numTol);
TP_actual = NaN(numSlices, numTol);

% Loop through the image slices
for slice=1:numSlices
    Idenoised = squeeze(MultiIdenoised(slice,:,:,:));
    Mask = squeeze(MultiMask(slice,:,:));
    Bgr = squeeze(MultiBgr(slice,:,:,:));
    IMasked = Idenoised .* Mask;
    
    xRes = size(Idenoised,1);
    zRes = size(Idenoised,2);
    numTimesteps = size(Idenoised,3);

    % Perform CaCLEAN on the simulated microscopy data
    clear CO;
    CO=CICRcleanSimp(Idenoised,Bgr,Mask,xyt_dim,'ApparentDiffusionK',60,'CleanDiffusionK',80,'CaCleanThreshold',10);
    % Get the CRU properties (segmented detected clusters)
    CruPropsCO = CRUProps(CO);
    numDetected = size(CruPropsCO.CRUProps,1);

    centroidsDetected = zeros(numDetected,2);
    clusterNumDetected = zeros(1,numDetected);
    detectedLabelMap = zeros(1,numDetected);
    uLabels = unique(CruPropsCO.CRULabel);
    % Get xy locations for ryr centroids
    for i=1:numDetected
        centroidsDetected(i,:) = CruPropsCO.CRUProps(i).WeightedCentroid;
        detectedLabelMap(i) = uLabels(i+1);
        clusterNumDetected(i) = i;
    end

    % Loop through the detection tolerances
    for tol=1:numTol
        availDetectedCentroids = centroidsDetected;
        tolClusterCenters = squeeze(MultiRyrClusterCenters(slice, tol, :, :));
        modeledClusterCenters = [tolClusterCenters(:,2), tolClusterCenters(:,1)]; % swapping z/x for normal imaging orientation
        numClusters = sum(~isnan(tolClusterCenters(:,1)));
        if numClusters > 0        
            realClusterAssocWithDetected = zeros(1,numDetected);
            truePositives = zeros(1,numClusters); % a 'hit'- correctly detected couplon
            falseNegatives = zeros(1,numClusters); % a 'miss'- couplon exists but not detected
            detectedCandidate = NaN(1,numClusters);
            detectedCandidateDist = NaN(1,numClusters);
            % First determine whether real (modeled) clusters fall within
            % which detected clusters
            for i=1:numClusters
                nearestPixel = [int16(modeledClusterCenters(i,2)),int16(modeledClusterCenters(i,1))];
                % Also look at surrounding pixels
                pixels = zeros(9,2);
                vals = zeros(9,1);
                n = 1;
                temp = [-1 0 1];
                for j=1:3
                    for k=1:3
                        pixels(n, 1) = temp(j) + nearestPixel(1);
                        pixels(n, 2) = temp(k) + nearestPixel(2);
                        if all(pixels(n, :) > 0)
                            vals(n,1) = CruPropsCO.CRULabel(pixels(n,1), pixels(n,2));
                        end
                        n = n+1;
                    end
                end
                uVals = unique(vals); 
                numUvals = size(uVals,1);
                dIdx = NaN(1,numUvals);
                dDist = NaN(1,numUvals);
                for j=1:numUvals
                    label = find(detectedLabelMap == uVals(j));
                    if ~isempty(label)
                        dIdx(j) = label;
                        dDist(j) = vecnorm(abs(centroidsDetected(dIdx(j),:) - modeledClusterCenters(i,:)),2,2);
                    end     
                end
                [dist, distIdx] = min(dDist);
                detectedIdx = dIdx(distIdx);
                if ~isnan(detectedIdx)
                    if ~isempty(detectedIdx)
                        detectedCandidate(i) = detectedIdx;
                        detectedCandidateDist(i) = vecnorm(abs(centroidsDetected(detectedIdx,:) - modeledClusterCenters(i,:)),2,2);
                    end
                end
            end
            % Adjust for duplicates
            i = 1;
            while i <= numClusters
               if isnan(detectedCandidate(1,i))
                   % false negative- no matching detected cluster
                   truePositives(i) = 0;
                   falseNegatives(i) = 1;
                   i = i + 1;
               else
                   % First get detected match for this real cluster
                   detectedIdx = detectedCandidate(i);
                   % Now also get a list of all other real clusters
                   % matching the same detected cluster
                   matches = find(detectedCandidate == detectedIdx);
                   numMatches = size(matches,2);
                   matchDist = detectedCandidateDist(matches);
                   % Now determine which real cluster is closer to the
                   % centroid of the matching detected cluster (the current
                   % real or an alternative)
                   [closestMatchDist, closestMatchIdx] = min(matchDist);
                   closestMatch = matches(closestMatchIdx);
                   if closestMatch == i
                       % The closest detected is also the closest match. Set as
                       % a true positive
                       truePositives(i) = detectedIdx;
                       falseNegatives(i) = 0;
                       % Remove from list of available (for plotting in fig2)
                       availDetectedCentroids(detectedIdx, :) = [0,0];                   
                       % Also remove this detected cluster from the other
                       % matches
                       for m = 1:numMatches
                           if m~= closestMatchIdx
                               detectedCandidateDist(matches(m)) = NaN;
                               detectedCandidate(matches(m)) = NaN;
                           end
                       end
                       i = i + 1;
                   else
                       % Another real cluster is closer to the current closest
                       % detected cluster. Remove this option from the list of
                       % options from the current real cluster list and repeat
                       % this iteration.
                       detectedCandidateDist(i) = NaN;
                       detectedCandidate(i) = NaN;
                   end
               end
            end

            TP = nnz(truePositives);
            FN = nnz(falseNegatives);
            FP = numDetected - TP;
            P = numClusters;

            Recall(slice, tol) = TP / P; % Sensitivity, recall, or TP rate
            Prec(slice, tol) = TP / (TP + FP); % Precision, positive predictive value
            FP_detected(slice, tol) = FP;
            FN_actual(slice, tol) = FN;
            TP_detected(slice, tol) = TP;
            TP_actual(slice, tol) = P;

            trues = zeros(TP,2);
            for i=1:numClusters
                if truePositives(i)
                    trues(i,:) = modeledClusterCenters(i,:);
                end
            end

            % set 0 values to nan for plotting
            trues(trues == 0) = NaN;
            availDetectedCentroids(availDetectedCentroids == 0) = NaN;

            % plots for single slice (example of detection/classification)
            if slice == fig3Slice
              if tol == fig3Tol
                fig3ModeledClusterCenters = modeledClusterCenters;
                fig3CO = CO;
                fig3CruPropsCruLabel = CruPropsCO.CRULabel;
                fig3Trues = trues;
                fig3CentroidsDetected = centroidsDetected;
                fig3AvailDetectedCentroids = availDetectedCentroids;
              end
            end
        end
        
    end
end


%% F I G U R E 3: Example of detection and classification
figure(1);
hold off;
%subplot(2,2,1)
colormap viridis;
%c = colorbar;
%c.Label.String = 'Fluorescence (a.u.)';
imagesc(squeeze(MultiIdenoised(fig3Slice, :, :, 4)));
set(gca,'YDir','normal');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
hold on
cbr = colorbar;
set(cbr,'YTick',0:2000:12000);
plot(fig3ModeledClusterCenters(:,1),fig3ModeledClusterCenters(:,2), 'm.','MarkerSize',10);
title('Fluorescence (a.u.)');

%hold off;
figure(2);
%subplot(2,2,2)
colormap viridis;
imagesc(fig3CO.CaRelease2D); caxis([0,50000]);
cbr = colorbar;
set(cbr,'YTick',0:10000:50000);
%imagesc(CO.CaRelease2D);
set(gca,'YDir','normal');
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
hold on;
plot(fig3ModeledClusterCenters(:,1),fig3ModeledClusterCenters(:,2), 'm.','MarkerSize',10);
title('CaCLEAN CRU map')

%hold off;
figure(3);
%subplot(2,2,3)
colormap gray;
imagesc(fig3CruPropsCruLabel); caxis([0,1])
set(gca,'YDir','normal')
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
title('Segmented detected clusters')
hold on
plot(fig3ModeledClusterCenters(:,1),fig3ModeledClusterCenters(:,2), 'm.','MarkerSize',10);

%hold off;
figure(4);
%subplot(2,2,4);
markerSize =  5;
%subplot(1,1,1)
hold on
plot(fig3ModeledClusterCenters(:,1),fig3ModeledClusterCenters(:,2), 'r.','MarkerSize',markerSize*3);
plot(fig3Trues(:,1),fig3Trues(:,2), 'g.','MarkerSize',markerSize*3);
plot(fig3CentroidsDetected(:,1),fig3CentroidsDetected(:,2),'bo','MarkerSize',markerSize*2,'LineWidth',3);
plot(fig3AvailDetectedCentroids(:,1),fig3AvailDetectedCentroids(:,2),'ro','MarkerSize',markerSize*2,'LineWidth',3);
axis([0 zRes+0.5 0 xRes+0.5]);
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
title('Classified detected clusters')
legend('FN (ground truth)', 'TP (ground truth)', 'TP (detected)', 'FP (detected)', 'Location','northwest');


%% F i g u r e  4:   Admissible window
figure(5);
hold on;
clf;
x = RyrTolerances * 1000;
fpPlot = shadedErrorBar(x,FP_detected,{@mean,@std},'lineprops','-r','patchSaturation',0.33);
fnPlot = shadedErrorBar(x,FN_actual,{@mean,@std},'lineprops','-b','patchSaturation',0.33);
tpdPlot = shadedErrorBar(x,TP_detected,{@mean,@std},'lineprops','-g','patchSaturation',0.33);
tpaPlot = shadedErrorBar(x,TP_actual,{@mean,@std},'lineprops','-m','patchSaturation',0.33);

tpaPlot.mainLine.LineWidth = 2;
tpdPlot.mainLine.LineWidth = 2;
fnPlot.mainLine.LineWidth = 2;
fpPlot.mainLine.LineWidth = 2;

grid on;
grid minor;
axis([0 1600 0 90]);
legend([tpaPlot.mainLine tpdPlot.mainLine fnPlot.mainLine fpPlot.mainLine], 'TP (ground truth)', 'TP (detected)', 'FN (ground truth)', 'FP (detected)', 'Location','northwest');
xlabel('Admissible window (nm)');
ylabel('Number of events');


%% F i g u r e   5:  Precision, recall, and F1-Score as a function of admissible window
figure(6);
clf;
hold on;
x = RyrTolerances * 1000;
snPlot = shadedErrorBar(x,Recall,{@mean,@std},'lineprops','-b','patchSaturation',0.33);
precPlot = shadedErrorBar(x,Prec,{@mean,@std},'lineprops','-r','patchSaturation',0.33);
snPlot.mainLine.LineWidth = 2;
precPlot.mainLine.LineWidth = 2;

meanPrec = squeeze(nanmean(Prec));
meanRecall = squeeze(nanmean(Recall));
f1Score = 2.0 .* (meanPrec .* meanRecall)./(meanPrec + meanRecall);
f1Plot = plot(x,f1Score,'-k','LineWidth',2);
[f1Max, f1MaxIdx] = max(f1Score);
maxF1X = x(f1MaxIdx);

grid on;
grid minor;
axis([0 1600 0 1.1]);
legend([snPlot.mainLine precPlot.mainLine f1Plot], 'Recall', 'Precision', 'F1-Score', 'Location','southeast');
title(plotTitle);
xlabel('Admissible window (nm)');

% Write figure to file
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
print(strcat(outfile, '.png'), '-dpng', '-r300');

% Save F1, Recall, and Precision values to text file 
%fullErrMinIdx = errMinIdx + startIdx - 1;
fileID = fopen(strcat(outfile,'.txt'),'w');
fprintf(fileID, strcat(plotTitle, '\n\n'));
outString = strcat('F1 Max admissible window:\t', string(maxF1X), 'nm\n');
fprintf(fileID, outString);
outString = strcat('F1 max:\t', string(f1Max), '\n');
fprintf(fileID, outString);
meanRecall = nanmean(Recall(:,f1MaxIdx));
stdRecall = std(Recall(:,f1MaxIdx));
outString = strcat('Recall :\t', string(meanRecall), '\t+/-\t', string(stdRecall) ,'\n');
fprintf(fileID, outString);
meanPrec = nanmean(Prec(:,f1MaxIdx));
stdPrec = std(Prec(:,f1MaxIdx));
outString = strcat('Precision :\t', string(meanPrec), '\t+/-\t', string(stdPrec) ,'\n\n');
fprintf(fileID, outString);

% Save intersection point F1, Recall and Precision values to text file
[diffPrecRecallIP,IPIdx] = min(abs(nanmean(Prec)-nanmean(Recall)));
outString = strcat('Intersection point admissible window:\t', string(x(IPIdx)), 'nm\n');
fprintf(fileID, outString);
outString = strcat('F1 score:\t', string(f1Score(IPIdx)), '\n');
fprintf(fileID, outString);
meanRecall = nanmean(Recall(:,IPIdx));
stdRecall = std(Recall(:,IPIdx));
outString = strcat('Recall :\t', string(meanRecall), '\t+/-\t', string(stdRecall) ,'\n');
fprintf(fileID, outString);
meanPrec = nanmean(Prec(:,IPIdx));
stdPrec = std(Prec(:,IPIdx));
outString = strcat('Precision :\t', string(meanPrec), '\t+/-\t', string(stdPrec) ,'\n\n');
fprintf(fileID, outString);

fclose(fileID);

%% F i g u r e   6:  Precision-recall plot
startNm = 80;
endNm = 1400;
nmInterval = 10;
startAW = startNm / nmInterval;
endAW = endNm / nmInterval;
numAW = round((endNm - startNm) / nmInterval);

x = squeeze(linspace(startNm, endNm - nmInterval, numAW))';
x2 = squeeze(linspace(startNm, endNm, numAW+1))';

dTPD = diff(TP_detected(:,startAW:endAW),1,2);
dTPGt = diff(TP_actual(:,startAW:endAW),1,2);

recall_all = nanmean(TP_detected(:,startAW:endAW) ./ TP_actual(:,startAW:endAW));
prec_all = nanmean(TP_detected(:,startAW:endAW) ./ (TP_detected(:,startAW:endAW) + FP_detected(:,startAW:endAW)));

figure(7);
plot(recall_all, prec_all, 'LineWidth',2);    
axis([0 1.0 0 1.1]);
xlabel('Recall');
ylabel('Precision');
title(plotTitle)
grid on;
grid minor;


%% F i g u r e   7:  Differential recall
ratio_dTPD = dTPD ./ dTPGt;
norm_dTPD = nanmean(dTPD ./ dTPGt);
sgf3_21 = sgolayfilt(norm_dTPD, 3, 21);
sgf_all = sgf3_21;    
f1 = fit(x, squeeze(sgf_all'),'gauss1');

figure(8);
hold on;
p1 = plot(f1,x,squeeze(sgf_all));
%set(p1, 'color', 'b');
xlabel('Distance z from focus (nm)');
yT = ylabel('Detection fraction: $\frac{\textrm{d\big(TP}_\textrm{detected}\big)}{\textrm{d\big(TP}_\textrm{ground truth}\big)}$','Interpreter','latex');
set(yT, 'FontSize', 12);
legend('Filtered results', 'Gaussian fit');
title('Fraction of detected sites as a function of distance z from focus');

%% Save workspace variables to file
%save(strcat(outfile,'.mat'));
