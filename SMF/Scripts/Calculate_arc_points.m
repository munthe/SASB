function [x,y] = Calculate_arc_points(ROC,ROCA,StartAngle,StopAngle,nr_sub_angle)

AngleStep = linspace(StartAngle,StopAngle,nr_sub_angle);


for k = 1:length(AngleStep)
%     arcx          = pi/2 + AngleStep(k);
    arcx          = AngleStep(k);
    fRad          = ROCA;
    x(k)          = fRad * cos(arcx);
    y(k)          = fRad * sin(arcx) - ROC;
end

if(sum(imag(x) ~= 0) ~= 0)
    error('Imaginary number in Calculate_arc_points')
end