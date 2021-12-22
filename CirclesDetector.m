tic
clear,clc
circleHoughTransform()
toc


function circleHoughTransform(minPoints, minRadius, maxRadius)

[image, m, n] = readImage('DottedCircles.png');
[labeledImage, labels] = labelImage();
setRadiiAndMinPoints(nargin);
radiiRange = maxRadius - minRadius + 1;
[theta, phase] = matchPhaseToEachRadius(radiiRange);
cosSinlinspace = cell(2, radiiRange);
accumulatorMat = zeros(m, n);
drawCirclesMat = ones(m, n);
accumulateOrDiminish(1, m, n)
while labels > 4
    [rad, xCen, yCen] = absAccumulator();
    relevantAreaMat = relevantArea(rad, xCen, yCen, minRadius, 10, m, n);
    [localAccumulatorMat, localMinRadius, radiiRangeLocal, localtheta, localX, localY] = localAccumulator(rad, minRadius, maxRadius, m, n);
    [~, xCen, yCen] = absAccumulator(1);
    [rad, xCen, yCen, maxVotes] = SetRadiusAndCenter(1, xCen, yCen, radiiRangeLocal, localMinRadius);
    if maxVotes < minPoints
        break;
    end
    relevantAreaMat = relevantArea(rad, xCen, yCen, minRadius, 2, m, n);
    accumulateOrDiminish(- 1, m, n)
    labeledImage(logical(relevantAreaMat)) = 1;
    [~, labels] = bwlabel(~ labeledImage, 4);
end
figure
drawCirclesMat = drawCirclesMat  - ~ labeledImage;
imshow(drawCirclesMat);


    function [image, m, n] = readImage(imageName)
        bwMat = image2BinaryMat(imageName);
        removeBordersAndGetSizes();

        function bwMat= image2BinaryMat(imageName)
            A = imread(imageName);
            A = im2gray(A);
            bwMat = imbinarize(A);
        end

        function removeBordersAndGetSizes()
            for side = 1:4
                while all(bwMat(:, 1))
                    bwMat(:, 1) = [];
                end
                bwMat = rot90(bwMat);
            end
            m = size(bwMat, 1);
            n = size(bwMat, 2);
            bwMat([1:2, m - 1:m], [1:2, n - 1:n]) = 1;
            image = bwMat;
        end
    end

    function [labeledImage, labels] = labelImage()
        CONNECTIVITY = 4;
        [image, labels] = bwlabel(~ image, CONNECTIVITY);
        meanPoints();

        function meanPoints()
            AVERAGE_HELPER = 0.3;
            LI = ones(m, n);
            for label = 1:labels
                [LIm, LIn] = find(image == label);
                LI(round(mean(LIm) - AVERAGE_HELPER):round(mean(LIm) + AVERAGE_HELPER), round(mean(LIn) - AVERAGE_HELPER):round(mean(LIn) + AVERAGE_HELPER)) = 0;
            end
            labeledImage = LI;
        end
    end

    function setRadiiAndMinPoints(numberOfInputArguments)
        switch numberOfInputArguments
            case 0
                minPoints = 3;
                RadiiDef();
            case 1
                RadiiDef();
        end
 
        function RadiiDef()
            minRadius = floor(m / 8) - 1;
            maxRadius = ceil(m / 4) + 1;
        end
    end

    function [theta, phase] = matchPhaseToEachRadius(radiiRange)
        theta = linspace(- pi + 2 * pi / radiiRange, pi, radiiRange);
        phase = complex(cos(theta),sin(theta));
    end

    function accumulateOrDiminish(sign, m, n)
        [y, x] = findBlackPixels(sign);
        for r = minRadius:maxRadius
            rCosSin = r - minRadius + 1;
            rphase = phase(rCosSin);
            if sign == 1
                fillCosSinlinspace(r, rCosSin)
            end
            createCircleAroundEachPixel(sign, r, rCosSin, rphase, m, n)
        end

        function [y, x] = findBlackPixels(sign)
            if sign == 1
                [y, x] = find(labeledImage == 0);
            else
                [y, x] = find(~ labeledImage .* relevantAreaMat == 1);
            end
        end

        function fillCosSinlinspace(r, rCosSin)
            omega = linspace(1 / r, 2 * pi, ceil(2 * pi * r));
            cosSinlinspace(:, rCosSin) = {cos(omega), sin(omega)};
        end

        function createCircleAroundEachPixel(sign, r, rCosSin, rphase, m, n)
            blackPoints = length(y);
            for blackPoint = 1:blackPoints
                a = round(x(blackPoint) - r .* cosSinlinspace{1, rCosSin});
                b = round(y(blackPoint) - r .* cosSinlinspace{2, rCosSin});
                keepInBordersPixelsOnly(m, n)
                addOrSubtractRadiusPhase(sign, rphase)
            end

            function keepInBordersPixelsOnly(m, n)
                log = a > 0 & a <= n & b > 0 & b <= m;
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

    function [rad, xCen, yCen] = absAccumulator(~)
        [A, AMax] = absAccumulatorMat(nargin);
        showAbsAccumulatorMat(0)
        findRadiusAndCenter(nargin, AMax)

        function [A, AMax] = absAccumulatorMat(numberOfInputArguments)
            switch numberOfInputArguments
                case 0
                    A = abs(accumulatorMat);
                otherwise
                    A = abs(localAccumulatorMat);
            end
            AMax = max(A(:));
        end

        function showAbsAccumulatorMat(binary)
            if binary
                figure
                imagesc(A)
                colormap hot
                axis equal off
            end
        end

        function findRadiusAndCenter(numberOfInputArguments, AMax)
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

    function relevantAreaMat = relevantArea(rad, xCen, yCen, minRadius, width, m, n)
        x = round(rad .* cosSinlinspace{1, rad - minRadius + 1} + xCen);
        y = round(rad .* cosSinlinspace{2, rad - minRadius + 1} + yCen);
        keepInBordersPixelsOnly(m, n)
        markRelevatAreaAndDrawCircle(length(x), width, m, n)

        function keepInBordersPixelsOnly(m, n)
            log = x > 0 & x <= n & y > 0 & y <= m;
            x = x(log);
            y = y(log);
        end

        function markRelevatAreaAndDrawCircle(lengthOfx, width, m, n)
            RA = zeros(m, n);
            for t = 1:lengthOfx
                xt = x(t);
                yt = y(t);
                RA(yt - min(yt - 1,width):yt + min(m - yt, width), xt - min(xt - 1, width):xt + min(n - xt, width)) = 1;
                if width == 2
                    drawCirclesMat(yt, xt) = 0;
                end
            end
            relevantAreaMat = RA;
        end
    end

    function [localAccumulatorMat, localMinRadius, localRadiiRange, localtheta, y, x] = localAccumulator(rad, minRadius, maxRadius, m, n)
        [y, x] = findBlackPixels();
        [localMinRadius, localMaxRadius, localRadiiRange] = setLocalRadii(rad, minRadius, maxRadius);
        [localtheta , localphase] = setLocalPhaseCoding(localRadiiRange);
        LAMat = zeros(m, n);
        for r = localMinRadius:localMaxRadius
            rCosSin = r - localMinRadius + 1;
            rPhase = localphase(rCosSin);
            createCircleAroundEachPixel(m, n, r, rCosSin, rPhase)
        end
        localAccumulatorMat = LAMat;

        function [y, x] = findBlackPixels()
            [y, x] = find(~ labeledImage .* relevantAreaMat == 1);
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

        function createCircleAroundEachPixel(m, n, r, rCosSin, rPhase)
            blackPoints = length(y);
            for blackPoint = 1:blackPoints
                a = round(x(blackPoint) - r .* cosSinlinspace{1, rCosSin});
                b = round(y(blackPoint) - r .* cosSinlinspace{2, rCosSin});
                keepInBordersPixelsOnly(m, n)
                addRadiusPhase(rPhase)
            end

            function keepInBordersPixelsOnly(m, n)
                log = a > 0 & a <= n & b > 0 & b <= m;
                a=a(log);
                b=b(log);
            end

            function addRadiusPhase(rPhase)
                circlePoints = length(a);
                for circlePoint = 1:circlePoints
                    bcp = b(circlePoint);
                    acp = a(circlePoint);
                    LAMat(bcp, acp) = LAMat(bcp, acp) + rPhase;
                end
            end
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