clear, clc, close all
tic
dottedCircles(2)
toc



function [radii, xCen, yCen, pointsSize] = setRadiiAndCenters( ...
        numOfcircles, scaleConst)
    radii = randi(round(scaleConst * [0.25, 0.5]), numOfcircles, 1);
    xCen = randi((scaleConst - 20) * [- 1, 1], numOfcircles, 1);
    yCen = randi((scaleConst - 20) * [- 1, 1], numOfcircles, 1);
    pointsSize = rand(numOfcircles, 1) * 2 + 1;
end

function [xNoise, yNoise, noiseSize] = setNoise(numOfNoisePoints, scaleConst)
    xNoise = randi((scaleConst - 2) * [- 1, 1], 1, numOfNoisePoints);
    yNoise = randi((scaleConst - 2) * [- 1, 1], 1, numOfNoisePoints);
    noiseSize = rand * 2 + 1;
end

function plotCorners(scaleConst)
    figure('visible', 'off')
    hold on
    plot(scaleConst * [- 1, 1, 1, - 1], scaleConst * [- 1, - 1, 1, 1], ...
        '.k', 'MarkerSize', 0.1)
    axis equal off
end

function outBordersPoints = rotateOutBordersPoints(rotMat, x, y, ...
        xCen, yCen, outBordersPoints, circle, SCALE_CONST)
    rotatedPoints = rotMat * [x(outBordersPoints) - xCen(circle); y(outBordersPoints) - yCen(circle)];
    x(outBordersPoints) = rotatedPoints(1, :) + xCen(circle);
    y(outBordersPoints) = rotatedPoints(2, :) + yCen(circle);
    outBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
end

function [x, y] = KeepInBordersPointsOnly(SCALE_CONST, minPoints, x, y, ...
            xCen, yCen, circle)
    outBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
    if nargin == 3
        while length(x) - length(outBordersPoints) < minPoints
            outBordersPoints = rotateOutBordersPoints([0, - 1; 1, 0], x, y, ...
            xCen, yCen, outBordersPoints, circle, SCALE_CONST);
        end
    end
    x(outBordersPoints) = [];
    y(outBordersPoints) = [];
end

function plotAndSaveCircles(whole ,imgName, pointsSize, radii, xCen, yCen, minPoints, numOfcircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
    plotCorners(SCALE_CONST)
    for circle = 1:numOfcircles
        if whole
            theta = linspace(1 / radii(circle), 2 * pi, 2 * pi * radii(circle));
        else
            theta = rand(1, randi([minPoints minPoints + 10])) * 2 * pi;
        end
        x = radii(circle) .* cos(theta) + xCen(circle);
        y = radii(circle) .* sin(theta) + yCen(circle);
        [x, y] = KeepInBordersPointsOnly(SCALE_CONST, minPoints, x, y, ...
            xCen, yCen, circle);
        plot(x, y, '.k', "MarkerSize", pointsSize(circle))
    end
    plot(xNoise, yNoise, '.k', "MarkerSize", noisePointsSize)
    saveas(gcf, imgName)
end

function mat = removeWhiteBorders(mat)
    for side = 1:4
        while all(mat(:, 1))
            mat(:, 1) = [];
        end
        mat = rot90(mat);
    end
    mat([1:2, end - 1:end], [1:2, end - 1:end]) = 1;
end

function imageShow(imageName)
    figure
    image = imbinarize(im2gray(imread([imageName,'.png'])));
    image = removeWhiteBorders(image);
    imshow(image)
end

function dottedCircles(numOfcircles, numOfNoisePoints, minPoints)
setInputArg(nargin);
SCALE_CONST = 200;
[radii, xCen, yCen, pointsSize] = setRadiiAndCenters(numOfcircles, SCALE_CONST);
[xNoise, yNoise, noisePointsSize] = setNoise(numOfNoisePoints, SCALE_CONST);
plotAndSaveCircles(false, 'DottedCircles.png', pointsSize, radii, xCen, yCen, minPoints, numOfcircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
plotAndSaveCircles(true, 'WholeCircles.png', pointsSize, radii, xCen, yCen, minPoints, numOfcircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
close all
imageShow('DottedCircles')
imageShow('WholeCircles')

    function setInputArg(numOfInputArguments)
        switch numOfInputArguments
            case 0
                numOfcircles = randi([3 10]);
                numOfNoisePoints = 0;
                minPoints = 30;
            case 1
                numOfNoisePoints = 0;
                minPoints = 30;
            case 2
                minPoints = max(30, round(3 * sqrt(numOfNoisePoints)));
        end
    end
end