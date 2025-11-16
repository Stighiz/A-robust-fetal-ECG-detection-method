function [] = plotData(x_axis, data, tit, x_label, y_label, trim)
    if nargin < 6
        trim = 0;
    end
    plot(x_axis, data);
    title(tit);
    xlabel(x_label);
    if trim==1
        xlim([0, 70]);
    end
    ylabel(y_label);
    grid on;
end