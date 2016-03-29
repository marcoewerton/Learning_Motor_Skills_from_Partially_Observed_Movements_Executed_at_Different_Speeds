function [] = getBoxSize()
% Return the *normalized* coordinates of the figure on the screen
% You should use this with 
%
% set(figure(1),'OuterPosition', b);
    
    try 
        set(gcf,'Units','normalized'); % this return the normalized value
        % set(gcf,'Units','pixels') 
        ii = get(gcf,'OuterPosition');
    end
    
    %fprintf('set(gcf,''Units'',''normalized''); set(gcf, ''OuterPosition'' , [%g %g %g %g]    );\n', ii);
    fprintf('Copy-paste:\n');
    fprintf('set_fig_position([%2.3g %2.3g %2.3g %2.3g]);\n', ii);
    
end