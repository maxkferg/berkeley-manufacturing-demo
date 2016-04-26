function showLoadPlot(XL,YL,ZL,SL)
% Show the load in X,Y,Z and S directions
    subplot (2,2,2);
    hold off
    plot(NaN)
    H = [XL; YL; ZL; SL];
    Nmax = numel(H);
    for i=1:Nmax
        h = bar(i, H(i));
        if i == 1, hold on, end
        if H(i) < 10
            col = 'g';
        elseif H(i) < 35
            col = 'y';
        else
            col = 'r';
        end
        set(h, 'FaceColor', col)
    end
    set(gca, 'XTickLabel', '')
    ylim([0 100])
    % Only perform these actions once
    title = get(get(gca, 'Title'), 'String');
    if isempty(title)
        setupLoadPlot(Nmax)
    end
end


function setupLoadPlot(Nmax)
    Xticks = {'X Load'; 'Y Load'; 'Z Load'; 'S Load'};
    xlabetxt = Xticks;
    ypos2 = -max(ylim)/50;
    text(1:Nmax,repmat(ypos2,Nmax,1),xlabetxt','horizontalalignment','center','verticalalignment','top','Rotation',0,'FontSize',6)
    title('Load as a % of the maximum rating')
    ylabel('% of maximum load capacity')
    set(gca,'FontSize',16)
end