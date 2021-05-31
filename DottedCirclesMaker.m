close all;clear,clc
DottedCircles()

function DottedCircles(circles,noise,MinPoints)

SetInputs(nargin);
SCALE_CONST = 200;
[Radii,Xcen,Ycen,PointsSize] = SetRadiiAndCenters(circles,SCALE_CONST);
[Xnoise,Ynoise,NoisePointsSize] = SetNoise(noise,SCALE_CONST);
PlotDottedCircels(circles,SCALE_CONST,Xnoise,Ynoise,NoisePointsSize)
PlotWholecircles(circles,SCALE_CONST,Xnoise,Ynoise,NoisePointsSize)
ImageShow('DottedCircles')
ImageShow('WholeCircles')


    function SetInputs(NumberOfInputArguments)
        switch NumberOfInputArguments
            case 0
                circles = randi([3 10]);
                noise = 0;
                MinPoints = 30;
            case 1
                noise = 0;
                MinPoints = 30;
            case 2
                MinPoints = max(30,round(3*sqrt(noise)));
        end
    end

    function [Radii,Xcen,Ycen,PointsSize] = SetRadiiAndCenters(circles,ScaleConst)
        Radii = randi(round(ScaleConst*[1/4,1/2]),circles,1);
        Xcen = randi((ScaleConst-20)*[-1,1],circles,1);
        Ycen = randi((ScaleConst-20)*[-1,1],circles,1);
        PointsSize = rand(circles,1)*2 + 1;
    end

    function [Xnoise,Ynoise,MarkerSizeNoise] = SetNoise(noise,ScaleConst)
        Xnoise = randi((ScaleConst - 2)*[-1,1],1,noise);
        Ynoise = randi((ScaleConst - 2)*[-1,1],1,noise);
        MarkerSizeNoise = rand*2 + 1;
    end

    function PlotDottedCircels(circles,ScaleConst,Xnoise,Ynoise,NoisePointsSize)
        PlotCorners(ScaleConst)
        for circle=1:circles
            t = rand(1,randi([MinPoints MinPoints + 10]))*2*pi;
            x = Radii(circle).*cos(t) + Xcen(circle);
            y = Radii(circle).*sin(t) + Ycen(circle);
            [x,y] = KeepInBordersPointsOnly(x,y,circle);
            plot(x,y,'.k',"MarkerSize",PointsSize(circle))
        end
        plot(Xnoise,Ynoise,'.k',"MarkerSize",NoisePointsSize)
        saveas(gcf,'DottedCircles.png')
    end

    function PlotWholecircles(circles,ScaleConst,Xnoise,Ynoise,NoisePointsSize)
        PlotCorners(ScaleConst)
        for circle = 1:circles
            t = linspace(1/Radii(circle),2*pi,2*pi*Radii(circle));
            x = Radii(circle).*cos(t) + Xcen(circle);
            y = Radii(circle).*sin(t) + Ycen(circle);
            [x,y] = KeepInBordersPointsOnly(x,y);
            plot(x,y,'.k',"MarkerSize",PointsSize(circle))
        end
        plot(Xnoise,Ynoise,'.k',"MarkerSize",NoisePointsSize)
        saveas(gcf,'WholeCircles.png')
        close all
    end

    function PlotCorners(ScaleConst)
        figure('visible','off')
        hold on
        plot(ScaleConst*[-1,1,1,-1],ScaleConst*[-1,-1,1,1],'.k','MarkerSize',0.1)
        axis equal off
    end

    function [x,y] = KeepInBordersPointsOnly(x,y,circle)
        OutBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
        if nargin==3
            while length(x) - length(OutBordersPoints) < MinPoints
                RotateOutBordersPoints([0,-1;1,0])
            end
        end
        
        x(OutBordersPoints) = [];
        y(OutBordersPoints) = [];
        function RotateOutBordersPoints(rot)
            z = rot*[x(OutBordersPoints) - Xcen(circle);y(OutBordersPoints) - Ycen(circle)];
            x(OutBordersPoints) = z(1,:) + Xcen(circle);
            y(OutBordersPoints) = z(2,:) + Ycen(circle);
            OutBordersPoints = find(abs(x) > SCALE_CONST|abs(y) > SCALE_CONST);
        end
    end

    function ImageShow(name)
        figure
        A = imread([name,'.png']);
        A = im2gray(A);
        A = imbinarize(A);
        A = RemoveBorders(A);
        imshow(A)
    end

    function rb = RemoveBorders(Mat)
        for side=1:4
            while all(Mat(:,1))
                Mat(:,1) = [];
            end
            Mat = rot90(Mat);
        end
        Mat([1:2,end-1:end],[1:2,end-1:end]) = 1;
        rb = Mat;
    end
end
