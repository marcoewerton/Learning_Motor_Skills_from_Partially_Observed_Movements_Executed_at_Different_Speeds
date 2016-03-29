clear variables;
close;
clc;

for i = 1:10
    
    close;
    figure(i);
    
    h = imfreehand('Closed', false);
    
    % get the position (x,y coordinates) of each point of the curve
    positions = getPosition(h);

end
