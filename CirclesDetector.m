tic
clear,clc
circleHoughTransform()
toc


function circleHoughTransform(minPoints, minRadius,maxRadius)

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
    [R, xCen, yCen] = absAccumulator();
    relevantAreaMat = relevantArea(R, xCen, yCen, minRadius, 10, m, n);
    [localAccumulatorMat, localMinRadius, radiiRangeLocal, localtheta, localX, localY] = localAccumulator(R, minRadius, maxRadius, m, n);
    [~, xCen, yCen] = absAccumulator(1);
    [R, xCen, yCen, MaxVotes] = SetRadiusAndCenter(1, xCen, yCen, radiiRangeLocal, localMinRadius);
    if MaxVotes < minPoints
        break;
    end
    relevantAreaMat = relevantArea(R, xCen, yCen, minRadius, 2, m, n);
    accumulateOrDiminish(-1, m, n)
    labeledImage(logical(relevantAreaMat)) = 1;
    [~, labels] = bwlabel(~ labeledImage, 4);
end
figure
drawCirclesMat = drawCirclesMat  - ~ labeledImage;
imshow(drawCirclesMat);


    function [image, m, n] = readImage(ImageName)
        bwMat = Image2BinaryMat(ImageName);
        RemoveBordersAndGetSizes();

        function bwMat= Image2BinaryMat(ImageName)
            A = imread(ImageName);
            A = im2gray(A);
            bwMat = imbinarize(A);
        end

        function RemoveBordersAndGetSizes()
            for side = 1:4
                while all(bwMat(:, 1))
                    bwMat(:, 1) = [];
                end
                bwMat = rot90(bwMat);
            end
            m = size(bwMat, 1);
            n = size(bwMat, 2);
            bwMat([1:2, m - 1:m],[1:2, n - 1:n]) = 1;
            image = bwMat;
        end
    end

    function [labeledImage, labels] = labelImage()
        CONNECTIVITY = 4;
        [image, labels] = bwlabel(~ image, CONNECTIVITY);
        MeanPoints();

        function MeanPoints()
            AVERAGE_HELPER = 0.3;
            LI = ones(m,n);
            for label = 1:labels
                [LIm, LIn] = find(image == label);
                LI(round(mean(LIm) - AVERAGE_HELPER):round(mean(LIm) + AVERAGE_HELPER), round(mean(LIn) - AVERAGE_HELPER):round(mean(LIn) + AVERAGE_HELPER)) = 0;
            end
            labeledImage = LI;
        end
    end

    function setRadiiAndMinPoints(NumberOfInputArguments)
        switch NumberOfInputArguments
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
        [y, x] = FindBlackPixels(sign);
        for r = minRadius:maxRadius
            rCosSin = r - minRadius + 1;
            rphase = phase(rCosSin);
            if sign == 1
                FillCosSinlinspace(r, rCosSin)
            end
            CreateCircleAroundEachPixel(sign, r, rCosSin, rphase, m, n)
        end

        function [y, x] = FindBlackPixels(sign)
            if sign == 1
                [y, x] = find(labeledImage == 0);
            else
                [y, x] = find(~ labeledImage .* relevantAreaMat == 1);
            end
        end

        function FillCosSinlinspace(r, rCosSin)
            omega = linspace(1 / r, 2 * pi, ceil(2 * pi * r));
            cosSinlinspace(:, rCosSin) = {cos(omega), sin(omega)};
        end

        function CreateCircleAroundEachPixel(sign, r, rCosSin, rphase, m, n)
            BlackPoints = length(y);
            for BlackPoint = 1:BlackPoints
                a = round(x(BlackPoint) - r .* cosSinlinspace{1, rCosSin});
                b = round(y(BlackPoint) - r .* cosSinlinspace{2, rCosSin});
                KeepInBordersPixelsOnly(m, n)
                AddOrSubtractRadiusPhase(sign, rphase)
            end

            function KeepInBordersPixelsOnly(m, n)
                log = a > 0 & a <= n & b > 0 & b <= m;
                a = a(log);
                b = b(log);
            end

            function AddOrSubtractRadiusPhase(sign, rphase)
                CirclePoints = length(a);
                for CirclePoint = 1:CirclePoints
                    bcp = b(CirclePoint);
                    acp = a(CirclePoint);
                    accumulatorMat(bcp, acp) = accumulatorMat(bcp, acp) + sign * rphase;
                end
            end
        end
    end

    function [R, xCen, yCen] = absAccumulator(~)
        [A, AMax] = absAccumulatorMat(nargin);
        ShowAbsAccumulatorMat(0)
        FindRadiusAndCenter(nargin, AMax)

        function [A, AMax] = absAccumulatorMat(NumberOfInputArguments)
            switch NumberOfInputArguments
                case 0
                    A = abs(accumulatorMat);
                otherwise
                    A = abs(localAccumulatorMat);
            end
            AMax = max(A(:));
        end

        function ShowAbsAccumulatorMat(binary)
            if binary
                figure
                imagesc(A)
                colormap hot
                axis equal off
            end
        end

        function FindRadiusAndCenter(NumberOfInputArguments, AMax)
            MIN_FRACT_OF_MAX_PIXEL = 0.8;
            [y, x] = MaxPixelVicinity(MIN_FRACT_OF_MAX_PIXEL, AMax);
            ExtractRadiusOutOfPhase(NumberOfInputArguments, length(y));
            xCen = round(mean(x));
            yCen = round(mean(y));

            function [y, x] = MaxPixelVicinity(MIN_FRACT_OF_MAX_PIXEL, AMax)
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

            function ExtractRadiusOutOfPhase(NumberOfInputArguments, NumberOfVicinityPixels)
                SumPhase = 0;
                if NumberOfInputArguments == 0
                    for VicinityPixel = 1:NumberOfVicinityPixels
                        SumPhase=SumPhase + accumulatorMat(y(VicinityPixel), x(VicinityPixel));
                    end
                    [~, PhaseRadius] = min(abs(theta - angle(SumPhase)));
                    R = PhaseRadius + minRadius - 1;
                else
                    for VicinityPixel = 1:NumberOfVicinityPixels
                        SumPhase = SumPhase + localAccumulatorMat(y(VicinityPixel), x(VicinityPixel));
                    end
                    [~, PhaseRadius] = min(abs(localtheta - angle(SumPhase)));
                    R = PhaseRadius + localMinRadius - 1;
                end
            end
        end
    end

    function relevantAreaMat = relevantArea(R, xCen, yCen, minRadius, width, m, n)
        x = round(R .* cosSinlinspace{1, R - minRadius + 1} + xCen);
        y = round(R .* cosSinlinspace{2, R - minRadius + 1} + yCen);
        KeepInBordersPixelsOnly(m, n)
        MarkRelevatAreaAndDrawCircle(length(x), width, m, n)

        function KeepInBordersPixelsOnly(m, n)
            log = x > 0 & x <= n & y > 0 & y <= m;
            x = x(log);
            y = y(log);
        end

        function MarkRelevatAreaAndDrawCircle(LengthOfx, width, m, n)
            RA = zeros(m, n);
            for t = 1:LengthOfx
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

    function [localAccumulatorMat, localMinRadius, LocalRadiiRange, localtheta, y, x] = localAccumulator(R, minRadius, maxRadius, m, n)
        [y, x] = FindBlackPixels();
        [localMinRadius, LocalMaxRadius, LocalRadiiRange] = SetLocalRadii(R, minRadius, maxRadius);
        [localtheta , Localphase] = SetLocalPhaseCoding(LocalRadiiRange);
        LAMat = zeros(m, n);
        for r = localMinRadius:LocalMaxRadius
            rCosSin = r - localMinRadius + 1;
            rPhase = Localphase(rCosSin);
            CreateCircleAroundEachPixel(m, n, r, rCosSin, rPhase)
        end
        localAccumulatorMat = LAMat;

        function [y, x] = FindBlackPixels()
            [y, x] = find(~ labeledImage .* relevantAreaMat == 1);
        end

        function [localMinRadius, LocalMaxRadius, LocalRadiiRange] = SetLocalRadii(R, minRadius, maxRadius)
            localMinRadius = max(R - 10, minRadius);
            LocalMaxRadius = min(R + 10, maxRadius);
            LocalRadiiRange = LocalMaxRadius - localMinRadius + 1;
        end

        function [localtheta, Localphase] = SetLocalPhaseCoding(LocalRadiiRange)
            localtheta = linspace(- pi + 2 * pi / LocalRadiiRange, pi, LocalRadiiRange);
            Localphase = complex(cos(localtheta), sin(localtheta));
        end

        function CreateCircleAroundEachPixel(m, n, r, rCosSin, rPhase)
            BlackPoints = length(y);
            for BlackPoint = 1:BlackPoints
                a = round(x(BlackPoint) - r .* cosSinlinspace{1, rCosSin});
                b = round(y(BlackPoint) - r .* cosSinlinspace{2, rCosSin});
                KeepInBordersPixelsOnly(m, n)
                AddRadiusPhase(rPhase)
            end

            function KeepInBordersPixelsOnly(m, n)
                log = a > 0 & a <= n & b > 0 & b <= m;
                a=a(log);
                b=b(log);
            end

            function AddRadiusPhase(rPhase)
                CirclePoints = length(a);
                for CirclePoint = 1:CirclePoints
                    bcp = b(CirclePoint);
                    acp = a(CirclePoint);
                    LAMat(bcp, acp) = LAMat(bcp, acp) + rPhase;
                end
            end
        end
    end

    function [Radius, xCenter, yCenter, MaxVotes] = SetRadiusAndCenter(STEP, xCen, yCen, radiiRangeLocal, localMinRadius)
        CenterAndRadius = zeros(2 * STEP + 1);
        BlackPoints = length(localX);
        for a = xCen - STEP:xCen + STEP
            for b = yCen - STEP:yCen + STEP
                RadiusAccumulator = zeros(1, radiiRangeLocal);
                for BlackPoint = 1:BlackPoints
                    AccumulateRadius(BlackPoint, a, b, localMinRadius, radiiRangeLocal);
                end
                AssignRadiusAndVotes(STEP, xCen, yCen)
            end
        end

        FindRadiusAndCenter(STEP, xCen, yCen, localMinRadius)

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
            CenterAndRadius(b - yCen + STEP + 1, a - xCen + STEP + 1) = complex(RA, maxRA);
        end

        function FindRadiusAndCenter(STEP, xCen, yCen, localMinRadius)
            RealCR = real(CenterAndRadius);
            ImagCR = imag(CenterAndRadius);
            MaxVotes = max(ImagCR(:));
            ExtractRadiusAndCenterOutOfCR(STEP, xCen, yCen, localMinRadius, MaxVotes)

            function ExtractRadiusAndCenterOutOfCR(STEP, xCen, yCen, localMinRadius, votes)
                RCR=round(mean(RealCR(ImagCR == votes)));
                [i, j] = find(ImagCR == votes);
                bCR = round(mean(i));
                aCR = round(mean(j));
                yCenter = yCen + bCR - STEP - 1;
                xCenter = xCen + aCR - STEP - 1;
                Radius = RCR + localMinRadius - 1;
            end
        end
    end
end