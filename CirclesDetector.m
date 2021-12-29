clear,clc,close all
tic
circleHoughTransform()
toc


function binaryMat = iamgeToBinaryMatrix(imgPath)
    binaryMat = imbinarize(im2gray(imread(imgPath)));
end

function mat = removeWhiteBorders(mat)
    for side = 1:4
        while all(mat(:, 1))
            mat(:, 1) = [];
        end
        mat = rot90(mat);
    end
end

function mat = whitePerimeter(mat, size1, size2)
    mat([1:2, size1 - 1:size1], [1:2, size2 - 1:size2]) = 1;
end

function [mat, size1, size2] = removeWhiteBordersAndGetSize(mat)
    mat = removeWhiteBorders(mat);
    [size1, size2] = size(mat);
    mat = whitePerimeter(mat, size1, size2);
end

function [mat, size1, size2] = readImageToMatrix(imgPath)
    binaryMat = iamgeToBinaryMatrix(imgPath);
    [mat, size1, size2] = removeWhiteBordersAndGetSize(binaryMat);
end

function labeledMat = centerLabelsPoints(size1, size2, labels, labeledMat)
    AVERAGE_HELPER = 0.3;
    tempMat = ones(size1, size2);
    for label = 1:labels
        [lebeledSpace1, labeledSpace2] = find(labeledMat == label);
        tempMat(round(mean(lebeledSpace1) - AVERAGE_HELPER):round(mean( ...
            lebeledSpace1) + AVERAGE_HELPER), round(mean(labeledSpace2) ...
            - AVERAGE_HELPER):round(mean(labeledSpace2) + ...
            AVERAGE_HELPER)) = 0;
    end
    labeledMat = tempMat;
end

function [labeledMat, labels] = connectedSpaces(mat, size1, size2)
    CONNECTIVITY = 4;
    [labeledMat, labels] = bwlabel(~ mat, CONNECTIVITY);
    labeledMat = centerLabelsPoints(size1, size2, labels, labeledMat);
end

function [minRad, maxRad] = radiiRangeDef(size1, size2)
    minRad = floor(min(size1, size2) / 8) - 1;
    maxRad = ceil(min(size1, size2) / 4) + 1;
end

function [theta, phase] = matchPhaseToEachRadius(radiiRange)
    theta = linspace(- pi + 2 * pi / radiiRange, pi, radiiRange);
    phase = complex(cos(theta),sin(theta));
end

function [y, x] = findRelevantBlackPixels(labeledMat, relAreaMat)
    switch nargin
        case 1
            [y, x] = find(labeledMat == 0);
        case 2
            [y, x] = find(~ labeledMat .* relAreaMat == 1);
    end
end

function cosSinlinspace = fillCosSinlinspace(r, rCosSin, cosSinlinspace)
    omega = linspace(1 / r, 2 * pi, ceil(2 * pi * r));
    cosSinlinspace(:, rCosSin) = {cos(omega), sin(omega)};
end

function [a, b] = keepInBordersPixelsOnly(a, b, size1, size2)
    logInsideBorders = a > 0 & a <= size2 & b > 0 & b <= size1;
    a = a(logInsideBorders);
    b = b(logInsideBorders);
end

function accMat = addOrSubtractRadiusPhase( ...
    sign, accMat, a, b, rPhase)
    crclPnts = length(a);
    for crclPnt = 1:crclPnts
        bCP = b(crclPnt);
        aCP = a(crclPnt);
        accMat(bCP, aCP) = accMat(bCP, aCP) + sign * rPhase;
    end
end

function accMat = createCircleAroundEachPixel(r, rCosSin, ...
    rPhase, size1, size2, x, y, cosSinlinspace, accMat, sign)
    blackPoints = length(y);
    for blackPoint = 1:blackPoints
        a = round(x(blackPoint) - r .* cosSinlinspace{1, rCosSin});
        b = round(y(blackPoint) - r .* cosSinlinspace{2, rCosSin});
        [a, b] = keepInBordersPixelsOnly(a, b, size1, size2);
        accMat = addOrSubtractRadiusPhase(sign, accMat, ...
            a, b, rPhase);
    end
end

function [cosSinlinspace, accMat] = fillAndCreate(minRad, ...
    maxRad, phase, cosSinlinspace, accMat, x, y, size1, size2, sign, fill)
    for r = minRad:maxRad
        rCosSin = r - minRad + 1;
        if fill
            cosSinlinspace = fillCosSinlinspace(r, rCosSin, cosSinlinspace);
        end
        rPhase = phase(rCosSin);
        accMat = createCircleAroundEachPixel(r, rCosSin, ...
            rPhase, size1, size2, x, y, cosSinlinspace, accMat, sign);
    end
end

function [cosSinlinspace, accMat] = accumulateOrDiminish(phase, ...
    size1, size2, labeledMat, cosSinlinspace, minRad, maxRad, ...
    accMat, sign, fill, relAreaMat)
    switch nargin
        case 10
            [y, x] = findRelevantBlackPixels(labeledMat);
        case 11
            [y, x] = findRelevantBlackPixels(labeledMat, relAreaMat);
    end    
    [cosSinlinspace, accMat] = fillAndCreate(minRad, ...
        maxRad, phase, cosSinlinspace, accMat, x, y, size1, ...
        size2, sign, fill);
end

function [absMat, maxAbsMat] = absoluteMatrix(mat)
    absMat = abs(mat);
    maxAbsMat = max(absMat(:));
end

function showAbsAccMat(binary, absMat)
    if binary
        figure
        imagesc(absMat)
        colormap hot
        axis equal off
    end
end

function [y, x] = keepOnlyVicinityPixels(y, x, yCen, xCen)
    VICINITY_RADIUS_SQUARED = 100;
    logInsideCircle = (y - yCen) .^ 2 + (x - xCen) .^ 2 < VICINITY_RADIUS_SQUARED;
    y = y(logInsideCircle);
    x = x(logInsideCircle);
end

function [y, x] = maxPixelVicinity(mat, MIN_FRACT_OF_MAX_PIXEL, matMax)
    [yMax, xMax] = find(mat == matMax);
    yMax = yMax(1);
    xMax = xMax(1);
    [y, x] = find(mat >= MIN_FRACT_OF_MAX_PIXEL * matMax);
    [y, x] = keepOnlyVicinityPixels(y, x, yMax, xMax);
end

function rad = extractRadiusOutOfPhase(numberOfVicinityPixels, ...
    accMat, x, y, theta, minRad)
    sumPhase = 0;
    for vicinityPixel = 1:numberOfVicinityPixels
        sumPhase = sumPhase + accMat(y(vicinityPixel), x(vicinityPixel));
    end
    [~, phaseRadius] = min(abs(theta - angle(sumPhase)));
    rad = phaseRadius + minRad - 1;
end

function [rad, xCen, yCen] = findRadAndCen(A, absAccMatMax, ...
    accMat, theta, minRad)
    MIN_FRACT_OF_MAX_PIXEL = 0.8;
    [y, x] = maxPixelVicinity(A, MIN_FRACT_OF_MAX_PIXEL, absAccMatMax);
    rad = extractRadiusOutOfPhase(length(y), accMat, x, y, theta, minRad);
    xCen = round(mean(x));
    yCen = round(mean(y));   
end

function [rad, xCen, yCen] = absAccumulator(accMat, theta, ...
    minRad)
    [absAccMat, absAccMatMax] = absoluteMatrix(accMat);
    showAbsAccMat(0, absAccMat)
    [rad, xCen, yCen] = findRadAndCen(absAccMat, absAccMatMax, accMat, ...
        theta, minRad);
end

function [x, y] = relevantPoints(rad, cosSinlinspace, xCen, yCen, minRad, size1, size2)
    x = round(rad .* cosSinlinspace{1, rad - minRad + 1} + xCen);
    y = round(rad .* cosSinlinspace{2, rad - minRad + 1} + yCen);
    [x, y] = keepInBordersPixelsOnly(x, y, size1, size2);
end

function relAreaMat =  markRelevatArea(x, y, xLength, width, ...
    size1, size2)
    tempMat = zeros(size1, size2);
    for t = 1:xLength
        x_t = x(t);
        y_t = y(t);
        tempMat(y_t - min(y_t - 1,width):y_t + min(size1 - y_t, width), ...
            x_t - min(x_t - 1, width):x_t + min(size2 - x_t, width)) = 1;
    end
    relAreaMat = tempMat;
end

function relAreaMat = relevantArea(rad, cosSinlinspace, xCen, yCen, minRad, width, size1, size2)
    [x, y] = relevantPoints(rad, cosSinlinspace, xCen, yCen, minRad, size1, size2);
    relAreaMat = markRelevatArea(x, y, length(x), width, size1, size2);
end

function [localminRad, localmaxRad, localRadiiRange] = setLocalRadii(rad, minRad, maxRad)
    RANGE_HELPER = 10;
    localminRad = max(rad - RANGE_HELPER, minRad);
    localmaxRad = min(rad + RANGE_HELPER, maxRad);
    localRadiiRange = localmaxRad - localminRad + 1;
end

function [localtheta, localphase] = setLocalPhaseCoding(localRadiiRange)
    localtheta = linspace(- pi + 2 * pi / localRadiiRange, pi, localRadiiRange);
    localphase = complex(cos(localtheta), sin(localtheta));
end

function [localAccumulatorMat, localminRad, localRadiiRange, localtheta, y, x] = localAccumulator( ...
        rad, minRad, maxRad, size1, size2, labeledMat, relAreaMat, cosSinlinspace)
    [y, x] = findRelevantBlackPixels(labeledMat, relAreaMat);
    [localminRad, localmaxRad, localRadiiRange] = setLocalRadii( ...
        rad, minRad, maxRad);
    [localtheta , localphase] = setLocalPhaseCoding(localRadiiRange);
    [~, localAccumulatorMat] = fillAndCreate(localminRad, ...
        localmaxRad, localphase, cosSinlinspace, zeros(size1, size2), x, y, size1, ...
        size2, 1, false);
end

function RadiusAccumulator = AccumulateRadius(bp, a, b,localminRad, radiiRangeLocal, RadiusAccumulator, localY, localX)
    ar = round(dist([localY(bp), localX(bp)], [a;b])) - localminRad + 1;
    if ar > 0 && ar <= radiiRangeLocal
        RadiusAccumulator(ar) = RadiusAccumulator(ar) + 1 ;
    end
end

function centerAndRadius = AssignRadiusAndVotes(STEP, xCen, yCen, a, b, RadiusAccumulator, centerAndRadius)
    maxRA = max(RadiusAccumulator);
    RA = find(RadiusAccumulator == maxRA);
    RA = round(mean(RA));
    centerAndRadius(b - yCen + STEP + 1, a - xCen + STEP + 1) = complex(RA, maxRA);
end

function [yCen, xCen, rad] = etractRadiusAndCenterOutOfCR(STEP, xCen, yCen, localminRad, maxVotes, realCR, imagCR)
    RCR = round(mean(realCR(imagCR == maxVotes)));
    [i, j] = find(imagCR == maxVotes);
    bCR = round(mean(i));
    aCR = round(mean(j));
    yCen = yCen + bCR - STEP - 1;
    xCen = xCen + aCR - STEP - 1;
    rad = RCR + localminRad - 1;
end

function [yCen, xCen, rad, maxVotes] = findRadiusAndCenter3(STEP, xCen, yCen, localminRad, centerAndRadius)
    realCR = real(centerAndRadius );
    imagCR = imag(centerAndRadius );
    maxVotes = max(imagCR(:));
    [yCen, xCen, rad] = etractRadiusAndCenterOutOfCR(STEP, xCen, yCen, localminRad, maxVotes, realCR, imagCR);
end

function centerAndRadius = findCenterAndRadius(xCen, yCen, STEP, radiiRangeLocal, blackPoints, localminRad, centerAndRadius, localY, localX)
    for a = xCen - STEP:xCen + STEP
        for b = yCen - STEP:yCen + STEP
            RadiusAccumulator = zeros(1, radiiRangeLocal);
            for blackPoint = 1:blackPoints
                RadiusAccumulator = AccumulateRadius(blackPoint, a, b, localminRad, radiiRangeLocal, RadiusAccumulator, localY, localX);
            end
            centerAndRadius = AssignRadiusAndVotes(STEP, xCen, yCen, a, b, RadiusAccumulator, centerAndRadius);
        end
    end
end

function [rad, xCen, yCen, maxVotes] = SetRadiusAndCenter(STEP, xCen, yCen, radiiRangeLocal, localminRad, localY, localX)
    centerAndRadius  = zeros(2 * STEP + 1);
    blackPoints = length(localX);
%     for a = xCen - STEP:xCen + STEP
%         for b = yCen - STEP:yCen + STEP
%             RadiusAccumulator = zeros(1, radiiRangeLocal);
%             for blackPoint = 1:blackPoints
%                 RadiusAccumulator = AccumulateRadius(blackPoint, a, b, localminRad, radiiRangeLocal, RadiusAccumulator, localY, localX);
%             end
%             centerAndRadius = AssignRadiusAndVotes(STEP, xCen, yCen, a, b, RadiusAccumulator, centerAndRadius);
%         end
%     end
    centerAndRadius = findCenterAndRadius(xCen, yCen, STEP, radiiRangeLocal, blackPoints, localminRad, centerAndRadius, localY, localX);
    [yCen, xCen, rad, maxVotes] = findRadiusAndCenter3(STEP, xCen, yCen, localminRad, centerAndRadius);
end

function [foundCirclesMat, relAreaMat] = markRelevatAreaAndDrawCircle(xLength, width, size1, size2, x, y, foundCirclesMat)
    RA = zeros(size1, size2);
    for t = 1:xLength
        xt = x(t);
        yt = y(t);
        RA(yt - min(yt - 1,width):yt + min(size1 - yt, width), xt - min(xt - 1, width):xt + min(size2 - xt, width)) = 1;
        if width == 2
            foundCirclesMat(yt, xt) = 0;
        end
    end
    relAreaMat = RA;
end

function [relAreaMat, foundCirclesMat] = relevantAreaAndDrawCircle(rad, xCen, yCen, minRad, width, size1, size2, cosSinlinspace, foundCirclesMat)
    [x, y] = relevantPoints(rad, cosSinlinspace, xCen, yCen, minRad, size1, size2);
    [foundCirclesMat, relAreaMat] = markRelevatAreaAndDrawCircle(length(x), width, size1, size2, x, y, foundCirclesMat);
end

function circleHoughTransform(minPoints, minRad, maxRad)
    [mat, size1, size2] = readImageToMatrix('DottedCircles.png');
    [labeledMat, labels] = connectedSpaces(mat, size1, size2);
    setRadiiAndMinPoints(nargin, size1, size2);
    radiiRange = maxRad - minRad + 1;
    [theta, phase] = matchPhaseToEachRadius(radiiRange);
    cosSinlinspace = cell(2, radiiRange);
    accMat = zeros(size1, size2);
    foundCirclesMat = ones(size1, size2);
    [cosSinlinspace, accMat] = accumulateOrDiminish(phase, size1, size2, ...
        labeledMat, cosSinlinspace, minRad, maxRad, accMat, 1, true);
    while labels > minPoints
        [rad, xCen, yCen] = absAccumulator(accMat, theta, minRad);
        relAreaMat = relevantArea(rad, cosSinlinspace, xCen, yCen, minRad, 10, ...
            size1, size2);
        [localAccumulatorMat, localminRad, radiiRangeLocal, localtheta, localX, localY] = localAccumulator( ...
            rad, minRad, maxRad, size1, size2, labeledMat, relAreaMat, cosSinlinspace);
        [~, xCen, yCen] = absAccumulator(localAccumulatorMat, localtheta, localminRad);
        [rad, xCen, yCen, maxVotes] = SetRadiusAndCenter(1, xCen, yCen, radiiRangeLocal, localminRad, localY, localX);
        if maxVotes < minPoints
            break;
        end
        [relAreaMat, foundCirclesMat]  = relevantAreaAndDrawCircle(rad, xCen, yCen, minRad, 2, size1, size2, cosSinlinspace, foundCirclesMat);
        
        [~, accMat] = accumulateOrDiminish(phase, size1, size2, ...
        labeledMat, cosSinlinspace, minRad, maxRad, accMat, -1, false, relAreaMat);
        labeledMat(logical(relAreaMat)) = 1;
        [~, labels] = bwlabel(~ labeledMat, 4);
    end
    figure
    foundCirclesMat = foundCirclesMat  - ~ labeledMat;
    imshow(foundCirclesMat);

        function setRadiiAndMinPoints(numberOfInputArguments, size1, size2)
        switch numberOfInputArguments
            case 0
                minPoints = 3;
                [minRad, maxRad] = radiiRangeDef(size1, size2);
            case 1
                [minRad, maxRad] = radiiRangeDef(size1, size2);
        end
    end
end