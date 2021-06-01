clear,clc
CircleHoughTransform()



function CircleHoughTransform(MinPoints,MinRadius,MaxRadius)

[image,m,n] = ReadImage('DottedCircles.png');
[LabeledImage,labels] = LabelImage();
SetRadiiAndMinPoints(nargin);
RadiiRange = MaxRadius - MinRadius + 1;
[theta,phase] = MatchPhaseToEachRadius(RadiiRange);
CosSinlinspace = cell(2,RadiiRange);
AccumulatorMat = zeros(m,n);
DrawCirclesMat = ones(m,n);
AccumulateOrDiminish(1,m,n)
while labels > 4
    [R,Xcen,Ycen] = AbsAccumulator();
    RelevantAreaMat = RelevantArea(R,Xcen,Ycen,MinRadius,10,m,n);
    [LocalAccumulatorMat,LocalMinRadius,RadiiRangeLocal,Localtheta,Localx,Localy] = LocalAccumulator(R,MinRadius,MaxRadius,m,n);
    [~,Xcen,Ycen] = AbsAccumulator(1);
    [R,Xcen,Ycen,MaxVotes] = SetRadiusAndCenter(1,Xcen,Ycen,RadiiRangeLocal,LocalMinRadius);
    if MaxVotes < MinPoints
        break;
    end
    RelevantAreaMat = RelevantArea(R,Xcen,Ycen,MinRadius,2,m,n);
    AccumulateOrDiminish(-1,m,n)
    LabeledImage(logical(RelevantAreaMat)) = 1;
    [~,labels] = bwlabel(~LabeledImage,4);
end
figure
DrawCirclesMat = DrawCirclesMat  - ~LabeledImage;
imshow(DrawCirclesMat);


    function [image,m,n] = ReadImage(ImageName)
        bwMat = Image2BinaryMat(ImageName);
        RemoveBordersAndGetSizes();

        function bwMat = Image2BinaryMat(ImageName)
            A = imread(ImageName);
            A = im2gray(A);
            bwMat = imbinarize(A);
        end

        function RemoveBordersAndGetSizes()
            for side = 1:4
                while all(bwMat(:,1))
                    bwMat(:,1) = [];
                end
                bwMat = rot90(bwMat);
            end
            m = size(bwMat,1);
            n = size(bwMat,2);
            bwMat([1:2,m - 1:m],[1:2,n - 1:n]) = 1;
            image = bwMat;
        end
    end

    function [LabeledImage,labels] = LabelImage()
        CONNECTIVITY = 4;
        [image,labels] = bwlabel(~image,CONNECTIVITY);
        MeanPoints();

        function MeanPoints()
            AVERAGE_HELPER = 0.3;
            LI = ones(m,n);
            for label = 1:labels
                [LIm,LIn] = find(image == label);
                LI(round(mean(LIm) - AVERAGE_HELPER):round(mean(LIm) + AVERAGE_HELPER),round(mean(LIn) - AVERAGE_HELPER):round(mean(LIn) + AVERAGE_HELPER)) = 0;
            end
            LabeledImage = LI;
        end
    end

    function SetRadiiAndMinPoints(NumberOfInputArguments)
        switch NumberOfInputArguments
            case 0
                MinPoints = 3;
                RadiiDef();
            case 1
                RadiiDef();
        end
 
        function RadiiDef()
            MinRadius = floor(m/8) - 1;
            MaxRadius = ceil(m/4) + 1;
        end
    end

    function [theta,phase] = MatchPhaseToEachRadius(RadiiRange)
        theta = linspace(-pi + 2*pi/RadiiRange,pi,RadiiRange);
        phase = complex(cos(theta),sin(theta));
    end

    function AccumulateOrDiminish(sign,m,n)
        [y,x] = FindBlackPixels(sign);
        for r = MinRadius:MaxRadius
            rCosSin = r - MinRadius + 1;
            rphase = phase(rCosSin);
            if sign == 1
                FillCosSinlinspace(r,rCosSin)
            end
            CreateCircleAroundEachPixel(sign,r,rCosSin,rphase,m,n)
        end

        function [y,x] = FindBlackPixels(sign)
            if sign == 1
                [y,x] = find(LabeledImage == 0);
            else
                [y,x] = find(~LabeledImage.*RelevantAreaMat == 1);
            end
        end

        function FillCosSinlinspace(r,rCosSin)
            omega = linspace(1/r,2*pi,ceil(2*pi*r));
            CosSinlinspace(:,rCosSin) = {cos(omega),sin(omega)};
        end

        function CreateCircleAroundEachPixel(sign,r,rCosSin,rphase,m,n)
            BlackPoints=length(y);
            for BlackPoint = 1:BlackPoints
                a = round(x(BlackPoint) - r.*CosSinlinspace{1,rCosSin});
                b = round(y(BlackPoint) - r.*CosSinlinspace{2,rCosSin});
                KeepInBordersPixelsOnly(m,n)
                AddOrSubtractRadiusPhase(sign,rphase)
            end

            function KeepInBordersPixelsOnly(m,n)
                log = a > 0 & a <= n & b > 0 & b <= m;
                a = a(log);
                b = b(log);
            end

            function AddOrSubtractRadiusPhase(sign,rphase)
                CirclePoints = length(a);
                for CirclePoint = 1:CirclePoints
                    bcp = b(CirclePoint);
                    acp = a(CirclePoint);
                    AccumulatorMat(bcp,acp) = AccumulatorMat(bcp,acp) + sign*rphase;
                end
            end
        end
    end

    function [R,Xcen,Ycen] = AbsAccumulator(~)
        [A,AMax] = AbsAccumulatorMat(nargin);
        ShowAbsAccumulatorMat(0)
        FindRadiusAndCenter(nargin,AMax)

        function [A,AMax] = AbsAccumulatorMat(NumberOfInputArguments)
            switch NumberOfInputArguments
                case 0
                    A = abs(AccumulatorMat);
                otherwise
                    A = abs(LocalAccumulatorMat);
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

        function FindRadiusAndCenter(NumberOfInputArguments,AMax)
            MIN_FRACT_OF_MAX_PIXEL = 0.8;
            [y,x] = MaxPixelVicinity(MIN_FRACT_OF_MAX_PIXEL,AMax);
            ExtractRadiusOutOfPhase(NumberOfInputArguments,length(y));
            Xcen = round(mean(x));
            Ycen = round(mean(y));

            function [y,x] = MaxPixelVicinity(MIN_FRACT_OF_MAX_PIXEL,AMax)
                [max_y,max_x] = find(A == AMax);
                max_y = max_y(1);
                max_x = max_x(1);
                [y,x] = find(A >= MIN_FRACT_OF_MAX_PIXEL*AMax);
                keepOnlyVicinityPixels()

                function keepOnlyVicinityPixels()
                    VICINITY_SQUARE_RADIUS = 100;
                    log = (y - max_y).^2 + (x - max_x).^2 < VICINITY_SQUARE_RADIUS;
                    y = y(log);
                    x = x(log);
                end
            end

            function ExtractRadiusOutOfPhase(NumberOfInputArguments,NumberOfVicinityPixels)
                SumPhase = 0;
                if NumberOfInputArguments == 0
                    for VicinityPixel = 1:NumberOfVicinityPixels
                        SumPhase = SumPhase + AccumulatorMat(y(VicinityPixel),x(VicinityPixel));
                    end
                    [~,PhaseRadius] = min(abs(theta - angle(SumPhase)));
                    R = PhaseRadius + MinRadius - 1;
                else
                    for VicinityPixel = 1:NumberOfVicinityPixels
                        SumPhase = SumPhase + LocalAccumulatorMat(y(VicinityPixel),x(VicinityPixel));
                    end
                    [~,PhaseRadius] = min(abs(Localtheta - angle(SumPhase)));
                    R = PhaseRadius + LocalMinRadius - 1;
                end
            end
        end
    end

    function RelevantAreaMat = RelevantArea(R,Xcen,Ycen,MinRadius,width,m,n)
        x = round(R.*CosSinlinspace{1,R - MinRadius + 1} + Xcen);
        y = round(R.*CosSinlinspace{2,R - MinRadius + 1} + Ycen);
        KeepInBordersPixelsOnly(m,n)
        MarkRelevatAreaAndDrawCircle(length(x),width,m,n)

        function KeepInBordersPixelsOnly(m,n)
            log = x > 0 & x <= n & y > 0 & y <= m;
            x = x(log);
            y = y(log);
        end

        function MarkRelevatAreaAndDrawCircle(LengthOfx,width,m,n)
            RA = zeros(m,n);
            for t = 1:LengthOfx
                xt = x(t);
                yt = y(t);
                RA(yt - min(yt - 1,width):yt + min(m - yt,width),xt - min(xt - 1,width):xt + min(n - xt,width)) = 1;
                if width == 2
                    DrawCirclesMat(yt,xt) = 0;
                end
            end
            RelevantAreaMat = RA;
        end
    end

    function [LocalAccumulatorMat,LocalMinRadius,LocalRadiiRange,Localtheta,y,x] = LocalAccumulator(R,MinRadius,MaxRadius,m,n)
        [y,x] = FindBlackPixels();
        [LocalMinRadius,LocalMaxRadius,LocalRadiiRange] = SetLocalRadii(R,MinRadius,MaxRadius);
        [Localtheta,Localphase] = SetLocalPhaseCoding(LocalRadiiRange);
        LAMat = zeros(m,n);
        for r = LocalMinRadius:LocalMaxRadius
            rCosSin = r - LocalMinRadius + 1;
            rPhase = Localphase(rCosSin);
            CreateCircleAroundEachPixel(m,n,r,rCosSin,rPhase)
        end
        LocalAccumulatorMat = LAMat;

        function [y,x] = FindBlackPixels()
            [y,x] = find(~LabeledImage.*RelevantAreaMat == 1);
        end

        function [LocalMinRadius,LocalMaxRadius,LocalRadiiRange] = SetLocalRadii(R,MinRadius,MaxRadius)
            LocalMinRadius = max(R - 10,MinRadius);
            LocalMaxRadius = min(R + 10,MaxRadius);
            LocalRadiiRange = LocalMaxRadius - LocalMinRadius + 1;
        end

        function [Localtheta,Localphase] = SetLocalPhaseCoding(LocalRadiiRange)
            Localtheta = linspace(-pi + 2*pi/LocalRadiiRange,pi,LocalRadiiRange);
            Localphase = complex(cos(Localtheta),sin(Localtheta));
        end

        function CreateCircleAroundEachPixel(m,n,r,rCosSin,rPhase)
            BlackPoints = length(y);
            for BlackPoint = 1:BlackPoints
                a = round(x(BlackPoint) - r.*CosSinlinspace{1,rCosSin});
                b = round(y(BlackPoint) - r.*CosSinlinspace{2,rCosSin});
                KeepInBordersPixelsOnly(m,n)
                AddRadiusPhase(rPhase)
            end

            function KeepInBordersPixelsOnly(m,n)
                log = a > 0 & a <= n & b > 0 & b <= m;
                a = a(log);
                b = b(log);
            end

            function AddRadiusPhase(rPhase)
                CirclePoints = length(a);
                for CirclePoint = 1:CirclePoints
                    bcp = b(CirclePoint);
                    acp = a(CirclePoint);
                    LAMat(bcp,acp) = LAMat(bcp,acp) + rPhase;
                end
            end
        end
    end

    function [Radius,Xcenter,Ycenter,MaxVotes] = SetRadiusAndCenter(STEP,Xcen,Ycen,RadiiRangeLocal,LocalMinRadius)
        CenterAndRadius = zeros(2*STEP + 1);
        BlackPoints = length(Localx);
        for a = Xcen - STEP:Xcen + STEP
            for b = Ycen - STEP:Ycen + STEP
                RadiusAccumulator = zeros(1,RadiiRangeLocal);
                for BlackPoint = 1:BlackPoints
                    AccumulateRadius(BlackPoint,a,b,LocalMinRadius,RadiiRangeLocal);
                end
                AssignRadiusAndVotes(STEP,Xcen,Ycen)
            end
        end

        FindRadiusAndCenter(STEP,Xcen,Ycen,LocalMinRadius)

        function AccumulateRadius(bp,a,b,LocalMinRadius,RadiiRangeLocal)
            ar = round(dist([Localy(bp),Localx(bp)],[a;b])) - LocalMinRadius + 1;
            if ar > 0 && ar <= RadiiRangeLocal
                RadiusAccumulator(ar) = RadiusAccumulator(ar) + 1 ;
            end
        end

       function AssignRadiusAndVotes(STEP,Xcen,Ycen)
            maxRA = max(RadiusAccumulator);
            RA = find(RadiusAccumulator == maxRA);
            RA = round(mean(RA));
            CenterAndRadius(b - Ycen + STEP + 1,a - Xcen + STEP + 1) = complex(RA,maxRA);
        end

        function FindRadiusAndCenter(STEP,Xcen,Ycen,LocalMinRadius)
            RealCR = real(CenterAndRadius);
            ImagCR = imag(CenterAndRadius);
            MaxVotes = max(ImagCR(:));
            ExtractRadiusAndCenterOutOfCR(STEP,Xcen,Ycen,LocalMinRadius,MaxVotes)

            function ExtractRadiusAndCenterOutOfCR(STEP,Xcen,Ycen,LocalMinRadius,votes)
                RCR = round(mean(RealCR(ImagCR == votes)));
                [i,j] = find(ImagCR == votes);
                bCR = round(mean(i));
                aCR = round(mean(j));
                Ycenter = Ycen + bCR - STEP - 1;
                Xcenter = Xcen + aCR - STEP - 1;
                Radius = RCR + LocalMinRadius - 1;
            end
        end
    end
end
