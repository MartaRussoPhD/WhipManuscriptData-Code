function FigEnlargeLabelForPng(pathname,fontsize)

f = gcf;
sp = gca;
sp = findobj('Type','Axes','Parent',f);

for isp = 1:length(sp)
    spp = sp(isp);
    
    spp.FontSize = fontsize;
    spp.XLabel.FontWeight = 'bold';
    spp.XLabel.FontSize = round(fontsize * 1.22);
    spp.YLabel.FontWeight = 'bold';
    spp.YLabel.FontSize = round(fontsize * 1.22);
    if iscell(spp.YLabel.String)
        for irow = 1:length(spp.YLabel.String)
            str = spp.YLabel.String{irow};
            maxstrlength = round(f.Position(4)*0.06);
            
            if length(str) > maxstrlength
                indspace = find(isspace(str));
                indendline = max(indspace(indspace <= maxstrlength)) - 1;
                newstr = str(indendline:end);
                spp.YLabel.String{irow} = str(1:indendline);
                spp.YLabel.String{irow,2} = newstr;
                
            end
        end
    else
        
        str = spp.YLabel.String;
        maxstrlength = round(f.Position(4)*0.06);
        if length(str) > maxstrlength
            indspace = find(isspace(str));
            indendline = max(indspace(indspace <= maxstrlength)) - 1;
            newstr = str(indendline:end);
            spp.YLabel.String = {str(1:indendline), newstr};
        end
        
    end
    
end
FigName = f.Name;

savefig(gcf,sprintf('%s/%s.fig',pathname,FigName));
saveas(gcf,sprintf('%s/%s.png',pathname,FigName));

close(gcf);

end