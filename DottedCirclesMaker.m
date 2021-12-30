clear, clc, close all
tic
dottedCircles(9)
toc



function [radii, xCen, yCen, pointsSize] = setRadiiAndCenters( ...
        circles, scaleConst)
    radii = randi(round(scaleConst * [0.25, 0.5]), circles, 1);
    xCen = randi((scaleConst - 20) * [- 1, 1], circles, 1);
    yCen = randi((scaleConst - 20) * [- 1, 1], circles, 1);
    pointsSize = rand(circles, 1) * 2 + 1;
end

function [xNoise,yNoise,MarkerSizeNoise] = setNoise(noise,scaleConst)
    xNoise = randi((scaleConst - 2)*[-1,1],1,noise);
    yNoise = randi((scaleConst - 2)*[-1,1],1,noise);
    MarkerSizeNoise = rand*2 + 1;
end








function dottedCircles(circles, noise, minPoints)

setInputs(nargin);
SCALE_CONST = 200;
[radii, xCen, yCen, pointsSize] = setRadiiAndCenters(circles, SCALE_CONST);
[xNoise, yNoise, noisePointsSize] = setNoise(noise, SCALE_CONST);
plotDottedCircels(circles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
plotWholecircles(circles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
imageShow('DottedCircles')
imageShow('WholeCircles')


    function setInputs(numberOfInputArguments)
        switch numberOfInputArguments
            case 0
                circles = randi([3 10]);
                noise = 0;
                minPoints = 30;
            case 1
                noise = 0;
                minPoints = 30;
            case 2
                minPoints = max(30, round(3 * sqrt(noise)));
        end
    end


    function plotDottedCircels(circles,scaleConst,xNoise,yNoise,noisePointsSize)
        PlotCorners(scaleConst)
        for circle=1:circles
            t = rand(1,randi([minPoints minPoints + 10]))*2*pi;
            x = radii(circle).*cos(t) + xCen(circle);
            y = radii(circle).*sin(t) + yCen(circle);
            [x,y] = KeepInBordersPointsOnly(x,y,circle);
            plot(x,y,'.k',"MarkerSize",pointsSize(circle))
        end
        plot(xNoise,yNoise,'.k',"MarkerSize",noisePointsSize)
        saveas(gcf,'DottedCircles.png')
    end

    function plotWholecircles(circles,scaleConst,xNoise,yNoise,noisePointsSize)
        PlotCorners(scaleConst)
        for circle = 1:circles
            t = linspace(1/radii(circle),2*pi,2*pi*radii(circle));
            x = radii(circle).*cos(t) + xCen(circle);
            y = radii(circle).*sin(t) + yCen(circle);
            [x,y] = KeepInBordersPointsOnly(x,y);
            plot(x,y,'.k',"MarkerSize",pointsSize(circle))
        end
        plot(xNoise,yNoise,'.k',"MarkerSize",noisePointsSize)
        saveas(gcf,'WholeCircles.png')
        close all
    end

    function PlotCorners(scaleConst)
        figure('visible','off')
        hold on
        plot(scaleConst*[-1,1,1,-1],scaleConst*[-1,-1,1,1],'.k','MarkerSize',0.1)
        axis equal off
    end

    function [x,y] = KeepInBordersPointsOnly(x,y,circle)
        OutBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
        if nargin==3
            while length(x) - length(OutBordersPoints) < minPoints
                RotateOutBordersPoints([0,-1;1,0])
            end
        end
        
        x(OutBordersPoints) = [];
        y(OutBordersPoints) = [];
        function RotateOutBordersPoints(rot)
            z = rot*[x(OutBordersPoints) - xCen(circle);y(OutBordersPoints) - yCen(circle)];
            x(OutBordersPoints) = z(1,:) + xCen(circle);
            y(OutBordersPoints) = z(2,:) + yCen(circle);
            OutBordersPoints = find(abs(x) > SCALE_CONST|abs(y) > SCALE_CONST);
        end
    end

    function imageShow(name)
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