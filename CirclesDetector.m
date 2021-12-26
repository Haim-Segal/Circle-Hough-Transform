clear,clc,close all
tic
circleHoughTransform()
toc



function binaryMat = imgToBinaryMat(imgPath)
    binaryMat = imbinarize(im2gray(imread(imgPath)));
end

function mat = removeBorders(mat)
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

function [mat, size1, size2] = removeBordersAndGetSizes(mat)
    mat = removeBorders(mat);
    [size1, size2] = size(mat);
    mat = whitePerimeter(mat, size1, size2);
end

function [mat, size1, size2] = readImage(imgPath)
    binaryMat = imgToBinaryMat(imgPath);
    [mat, size1, size2] = removeBordersAndGetSizes(binaryMat);
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

function [theta, phase] = matchPhaseToEachRadius(radiiRange)
    theta = linspace(- pi + 2 * pi / radiiRange, pi, radiiRange);
    phase = complex(cos(theta),sin(theta));
end

function [y, x] = findBlackPixels1(labeledMat)
    [y, x] = find(labeledMat == 0);
end

function cosSinlinspace = fillCosSinlinspace(r, rCosSin, cosSinlinspace)
    omega = linspace(1 / r, 2 * pi, ceil(2 * pi * r));
    cosSinlinspace(:, rCosSin) = {cos(omega), sin(omega)};
end

function [a, b] = keepInBordersPixelsOnly1(a, b, size1, size2)
    log = a > 0 & a <= size2 & b > 0 & b <= size1;
    a = a(log);
    b = b(log);
end

function accumulatorMat = addOrSubtractRadiusPhase( ...
        accumulatorMat, a, b, rphase)
    circlePoints = length(a);
    for circlePoint = 1:circlePoints
        bcp = b(circlePoint);
        acp = a(circlePoint);
        accumulatorMat(bcp, acp) = accumulatorMat( ...
            bcp, acp) + rphase;
    end
end

function accumulatorMat = createCircleAroundEachPixel1(accumulatorMat, ...
    x, y, r, rCosSin, cosSinlinspace, rphase, size1, size2)
        blackPoints = length(y);
        for blackPoint = 1:blackPoints
            a = round(x(blackPoint) - r .* cosSinlinspace{1, rCosSin});
            b = round(y(blackPoint) - r .* cosSinlinspace{2, rCosSin});
            [a, b] = keepInBordersPixelsOnly1(a, b, size1, size2);
            accumulatorMat = addOrSubtractRadiusPhase(accumulatorMat, ...
                a, b, rphase);
        end
end

function [cosSinlinspace, accumulatorMat] = fillAndCreate(minRadius, ...
    maxRadius, phase, cosSinlinspace, accumulatorMat, x, y, size1, size2)
    for r = minRadius:maxRadius
        rCosSin = r - minRadius + 1;
        rphase = phase(rCosSin);
        cosSinlinspace = fillCosSinlinspace(r, rCosSin, cosSinlinspace);
        accumulatorMat = createCircleAroundEachPixel1(accumulatorMat, ...
            x, y, r, rCosSin, cosSinlinspace, rphase, size1, size2);
    end
end

function [cosSinlinspace, accumulatorMat] = accumulate(phase, size1, ...
    size2, labeledMat, cosSinlinspace, minRadius, maxRadius, ...
    accumulatorMat)
    [y, x] = findBlackPixels1(labeledMat);
    [cosSinlinspace, accumulatorMat] = fillAndCreate(minRadius, ...
        maxRadius, phase, cosSinlinspace,accumulatorMat, x, y, size1, ...
        size2);
end

function [A, AMax] = absAccumulatorMat1(accumulatorMat)
    A = abs(accumulatorMat);
    AMax = max(A(:));
end

function showAbsAccumulatorMat(binary, A)
    if binary
        figure
        imagesc(A)
        colormap hot
        axis equal off
    end
end

function [y, x] = keepOnlyVicinityPixels1(y, x, max_y, max_x)
    VICINITY_SQUARE_RADIUS = 100;
    log = (y - max_y) .^ 2 + (x - max_x) .^ 2 < VICINITY_SQUARE_RADIUS;
    y = y(log);
    x = x(log);
end

function [y, x] = maxPixelVicinity1(A, MIN_FRACT_OF_MAX_PIXEL, AMax)
    [max_y, max_x] = find(A == AMax);
    max_y = max_y(1);
    max_x = max_x(1);
    [y, x] = find(A >= MIN_FRACT_OF_MAX_PIXEL * AMax);
    [y, x] = keepOnlyVicinityPixels1(y, x, max_y, max_x);
end

function rad = extractRadiusOutOfPhase1(accumulatorMat, x, y, ...
    numberOfVicinityPixels, theta, minRadius)
    sumPhase = 0;
    for vicinityPixel = 1:numberOfVicinityPixels
        sumPhase=sumPhase + accumulatorMat(y(vicinityPixel), x(vicinityPixel));
    end
    [~, phaseRadius] = min(abs(theta - angle(sumPhase)));
    rad = phaseRadius + minRadius - 1;
end

function [rad, xCen, yCen] = findRadiusAndCenter1(A, AMax, ...
            accumulatorMat, theta, minRadius)
    MIN_FRACT_OF_MAX_PIXEL = 0.8;
    [y, x] = maxPixelVicinity1(A, MIN_FRACT_OF_MAX_PIXEL, AMax);
    rad = extractRadiusOutOfPhase1(accumulatorMat, x, y, length(y), theta, minRadius);
    xCen = round(mean(x));
    yCen = round(mean(y));   
end

function [rad, xCen, yCen] = absAccumulator1(accumulatorMat, theta, ...
    minRadius)
    [A, AMax] = absAccumulatorMat1(accumulatorMat);
    showAbsAccumulatorMat(0, A)
    [rad, xCen, yCen] = findRadiusAndCenter1(A, AMax, accumulatorMat, ...
        theta, minRadius);
end

function [x, y] = keepInBordersPixelsOnly2(x, y, size1, size2)
    log = x > 0 & x <= size2 & y > 0 & y <= size1;
    x = x(log);
    y = y(log);
end

function relevantAreaMat =  markRelevatArea(x, y, xLength, width, ...
    size1, size2)
    tempMat = zeros(size1, size2);
    for t = 1:xLength
        x_t = x(t);
        y_t = y(t);
        tempMat(y_t - min(y_t - 1,width):y_t + min(size1 - y_t, width), ...
            x_t - min(x_t - 1, width):x_t + min(size2 - x_t, width)) = 1;
    end
    relevantAreaMat = tempMat;
end

function relevantAreaMat = relevantArea1(rad, cosSinlinspace, xCen, yCen, minRadius, width, size1, size2)
    x = round(rad .* cosSinlinspace{1, rad - minRadius + 1} + xCen);
    y = round(rad .* cosSinlinspace{2, rad - minRadius + 1} + yCen);
    [x, y] = keepInBordersPixelsOnly2(x, y, size1, size2);
    relevantAreaMat = markRelevatArea(x, y, length(x), width, size1, size2);
end

function [y, x] = findBlackPixels3(labeledMat, relevantAreaMat)
    [y, x] = find(~ labeledMat .* relevantAreaMat == 1);
end

function [localMinRadius, localMaxRadius, localRadiiRange] = setLocalRadii(rad, minRadius, maxRadius)
    localMinRadius = max(rad - 10, minRadius);
    localMaxRadius = min(rad + 10, maxRadius);
    localRadiiRange = localMaxRadius - localMinRadius + 1;
end

function [localtheta, localphase] = setLocalPhaseCoding(localRadiiRange)
    localtheta = linspace(- pi + 2 * pi / localRadiiRange, pi, localRadiiRange);
    localphase = complex(cos(localtheta), sin(localtheta));
end

function [a, b] = keepInBordersPixelsOnly5(a, b, size1, size2)
    log = a > 0 & a <= size2 & b > 0 & b <= size1;
    a=a(log);
    b=b(log);
end

function LAMat = addRadiusPhase(a, b, rPhase, LAMat)
    circlePoints = length(a);
    for circlePoint = 1:circlePoints
        bcp = b(circlePoint);
        acp = a(circlePoint);
        LAMat(bcp, acp) = LAMat(bcp, acp) + rPhase;
    end
end

function LAMat = createCircleAroundEachPixel2(x, y, cosSinlinspace, ...
        size1, size2, r, rCosSin, rPhase, LAMat)
       
    blackPoints = length(y);
    for blackPoint = 1:blackPoints
        a = round(x(blackPoint) - r .* cosSinlinspace{1, rCosSin});
        b = round(y(blackPoint) - r .* cosSinlinspace{2, rCosSin});
        [a, b] = keepInBordersPixelsOnly5(a, b, size1, size2);
        LAMat = addRadiusPhase(a, b, rPhase, LAMat);
    end
end

function [localAccumulatorMat, localMinRadius, localRadiiRange, localtheta, y, x] = localAccumulator( ...
        rad, minRadius, maxRadius, size1, size2, labeledMat, relevantAreaMat, cosSinlinspace)
    [y, x] = findBlackPixels3(labeledMat, relevantAreaMat);
    [localMinRadius, localMaxRadius, localRadiiRange] = setLocalRadii( ...
        rad, minRadius, maxRadius);
    [localtheta , localphase] = setLocalPhaseCoding(localRadiiRange);
    LAMat = zeros(size1, size2);
    for r = localMinRadius:localMaxRadius
        rCosSin = r - localMinRadius + 1;
        rPhase = localphase(rCosSin);
        LAMat = createCircleAroundEachPixel2(x, y, cosSinlinspace, size1, ...
            size2, r, rCosSin, rPhase, LAMat);
    end
    localAccumulatorMat = LAMat;
end




function circleHoughTransform(minPoints, minRadius, maxRadius)

[mat, size1, size2] = readImage('DottedCircles.png');
[labeledMat, labels] = connectedSpaces(mat, size1, size2);
setRadiiAndMinPoints(nargin);
radiiRange = maxRadius - minRadius + 1;
[theta, phase] = matchPhaseToEachRadius(radiiRange);
cosSinlinspace = cell(2, radiiRange);
accumulatorMat = zeros(size1, size2);
drawCirclesMat = ones(size1, size2);
[cosSinlinspace, accumulatorMat] = accumulate(phase, size1, size2, ...
    labeledMat, cosSinlinspace, minRadius, maxRadius, accumulatorMat);


while labels > 4
    [rad, xCen, yCen] = absAccumulator1(accumulatorMat, theta, minRadius);
    relevantAreaMat = relevantArea1(rad, cosSinlinspace, xCen, yCen, minRadius, 10, ...
        size1, size2);
    [localAccumulatorMat, localMinRadius, radiiRangeLocal, localtheta, localX, localY] = localAccumulator( ...
        rad, minRadius, maxRadius, size1, size2, labeledMat, relevantAreaMat, cosSinlinspace);
    
    [~, xCen, yCen] = absAccumulator2(localAccumulatorMat);
    
    [rad, xCen, yCen, maxVotes] = SetRadiusAndCenter(1, xCen, yCen, radiiRangeLocal, localMinRadius);
    if maxVotes < minPoints
        break;
    end
    relevantAreaMat = relevantArea(rad, xCen, yCen, minRadius, 2, size1, size2);
    diminish(- 1, size1, size2, labeledMat)
    labeledMat(logical(relevantAreaMat)) = 1;
    [~, labels] = bwlabel(~ labeledMat, 4);
end
figure
drawCirclesMat = drawCirclesMat  - ~ labeledMat;
imshow(drawCirclesMat);

    function setRadiiAndMinPoints(numberOfInputArguments)
        switch numberOfInputArguments
            case 0
                minPoints = 3;
                RadiiDef();
            case 1
                RadiiDef();
        end
 
        function RadiiDef()
            minRadius = floor(size1 / 8) - 1;
            maxRadius = ceil(size1 / 4) + 1;
        end
    end

  

    function diminish(sign, size1, size2, labeledMat)
        [y, x] = findBlackPixels2(sign, labeledMat);
        for r = minRadius:maxRadius
            rCosSin = r - minRadius + 1;
            rphase = phase(rCosSin);
            if sign == 1
                fillCosSinlinspace(r, rCosSin)
            end
            createCircleAroundEachPixel(sign, r, rCosSin, rphase, size1, size2)
        end

        function [y, x] = findBlackPixels2(sign, labeledMat)
            if sign == 1
                [y, x] = find(labeledMat == 0);
            else
                [y, x] = find(~ labeledMat .* relevantAreaMat == 1);
            end
        end

        function fillCosSinlinspace(r, rCosSin)
            omega = linspace(1 / r, 2 * pi, ceil(2 * pi * r));
            cosSinlinspace(:, rCosSin) = {cos(omega), sin(omega)};
        end

        function createCircleAroundEachPixel(sign, r, rCosSin, rphase, size1, size2)
            blackPoints = length(y);
            for blackPoint = 1:blackPoints
                a = round(x(blackPoint) - r .* cosSinlinspace{1, rCosSin});
                b = round(y(blackPoint) - r .* cosSinlinspace{2, rCosSin});
                keepInBordersPixelsOnly3(size1, size2)
                addOrSubtractRadiusPhase(sign, rphase)
            end

            function keepInBordersPixelsOnly3(size1, size2)
                log = a > 0 & a <= size2 & b > 0 & b <= size1;
                a = a(log);
                b = b(log);
            end

            function addOrSubtractRadiusPhase(sign, rphase)
                circlePoints = length(a);
                for circlePoint = 1:circlePoints
                    bcp = b(circlePoint);
                    acp = a(circlePoint);
                    accumulatorMat(bcp, acp) = accumulatorMat(bcp, acp) + sign * rphase;
                end
            end
        end
    end

    function [rad, xCen, yCen] = absAccumulator2(localAccumulatorMat)
        [A, AMax] = absAccumulatorMat2(localAccumulatorMat);
        showAbsAccumulatorMat(0, A)
        findRadiusAndCenter1(nargin, AMax)

        function [A, AMax] = absAccumulatorMat2(localAccumulatorMat)
            A = abs(localAccumulatorMat);  
            AMax = max(A(:));
        end

        function findRadiusAndCenter1(numberOfInputArguments, AMax)
            MIN_FRACT_OF_MAX_PIXEL = 0.8;
            [y, x] = maxPixelVicinity(MIN_FRACT_OF_MAX_PIXEL, AMax);
            extractRadiusOutOfPhase(numberOfInputArguments, length(y));
            xCen = round(mean(x));
            yCen = round(mean(y));

            function [y, x] = maxPixelVicinity(MIN_FRACT_OF_MAX_PIXEL, AMax)
                [max_y, max_x] = find(A == AMax);
                max_y = max_y(1);
                max_x = max_x(1);
                [y, x] = find(A >= MIN_FRACT_OF_MAX_PIXEL * AMax);
                keepOnlyVicinityPixels()

                function keepOnlyVicinityPixels()
                    VICINITY_SQUARE_RADIUS = 100;
                    log = (y - max_y) .^ 2 + (x - max_x) .^ 2 < VICINITY_SQUARE_RADIUS;
                    y = y(log);
                    x = x(log);
                end
            end

            function extractRadiusOutOfPhase(numberOfInputArguments, numberOfVicinityPixels)
                sumPhase = 0;
                if numberOfInputArguments == 0
                    for vicinityPixel = 1:numberOfVicinityPixels
                        sumPhase=sumPhase + accumulatorMat(y(vicinityPixel), x(vicinityPixel));
                    end
                    [~, phaseRadius] = min(abs(theta - angle(sumPhase)));
                    rad = phaseRadius + minRadius - 1;
                else
                    for vicinityPixel = 1:numberOfVicinityPixels
                        sumPhase = sumPhase + localAccumulatorMat(y(vicinityPixel), x(vicinityPixel));
                    end
                    [~, phaseRadius] = min(abs(localtheta - angle(sumPhase)));
                    rad = phaseRadius + localMinRadius - 1;
                end
            end
        end
    end

    function relevantAreaMat = relevantArea(rad, xCen, yCen, minRadius, width, size1, size2)
        x = round(rad .* cosSinlinspace{1, rad - minRadius + 1} + xCen);
        y = round(rad .* cosSinlinspace{2, rad - minRadius + 1} + yCen);
        keepInBordersPixelsOnly4(size1, size2)
        markRelevatAreaAndDrawCircle(length(x), width, size1, size2)

        function keepInBordersPixelsOnly4(size1, size2)
            log = x > 0 & x <= size2 & y > 0 & y <= size1;
            x = x(log);
            y = y(log);
        end

        function markRelevatAreaAndDrawCircle(xLength, width, size1, size2)
            RA = zeros(size1, size2);
            for t = 1:xLength
                xt = x(t);
                yt = y(t);
                RA(yt - min(yt - 1,width):yt + min(size1 - yt, width), xt - min(xt - 1, width):xt + min(size2 - xt, width)) = 1;
                if width == 2
                    drawCirclesMat(yt, xt) = 0;
                end
            end
            relevantAreaMat = RA;
        end
    end





    


    function [radius, xCenter, yCenter, maxVotes] = SetRadiusAndCenter(STEP, xCen, yCen, radiiRangeLocal, localMinRadius)
        centerAndRadius  = zeros(2 * STEP + 1);
        blackPoints = length(localX);
        for a = xCen - STEP:xCen + STEP
            for b = yCen - STEP:yCen + STEP
                RadiusAccumulator = zeros(1, radiiRangeLocal);
                for blackPoint = 1:blackPoints
                    AccumulateRadius(blackPoint, a, b, localMinRadius, radiiRangeLocal);
                end
                AssignRadiusAndVotes(STEP, xCen, yCen)
            end
        end

        findRadiusAndCenter(STEP, xCen, yCen, localMinRadius)

        function AccumulateRadius(bp, a, b,localMinRadius,radiiRangeLocal)
            ar = round(dist([localY(bp), localX(bp)], [a;b])) - localMinRadius + 1;
            if ar > 0 && ar <= radiiRangeLocal
                RadiusAccumulator(ar) = RadiusAccumulator(ar) + 1 ;
            end
        end

       function AssignRadiusAndVotes(STEP, xCen, yCen)
            maxRA = max(RadiusAccumulator);
            RA = find(RadiusAccumulator == maxRA);
            RA = round(mean(RA));
            centerAndRadius (b - yCen + STEP + 1, a - xCen + STEP + 1) = complex(RA, maxRA);
        end

        function findRadiusAndCenter(STEP, xCen, yCen, localMinRadius)
            realCR = real(centerAndRadius );
            imagCR = imag(centerAndRadius );
            maxVotes = max(imagCR(:));
            etractRadiusAndCenterOutOfCR(STEP, xCen, yCen, localMinRadius, maxVotes)

            function etractRadiusAndCenterOutOfCR(STEP, xCen, yCen, localMinRadius, votes)
                RCR = round(mean(realCR(imagCR == votes)));
                [i, j] = find(imagCR == votes);
                bCR = round(mean(i));
                aCR = round(mean(j));
                yCenter = yCen + bCR - STEP - 1;
                xCenter = xCen + aCR - STEP - 1;
                radius = RCR + localMinRadius - 1;
            end
        end
    end
end