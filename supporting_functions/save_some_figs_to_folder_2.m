function [] = save_some_figs_to_folder_2(save_folder, save_name, fig_vec, fig_type)

if not(exist(save_folder,'dir'))
    mkdir(save_folder)
else
    warning('folder already exists')
end

if isempty(fig_vec)
    figHandles = findobj('Type', 'figure')
    for i_f = 1:length(figHandles)
        fig_vec(i_f) = figHandles(i_f).Number;
    end
end

if isempty(fig_type)
    fig_type = {'fig','svg','png'};
end
    
h = get(0,'children');
h = flipud(h);

for i=fig_vec
    set(i,'PaperPositionMode','auto')  
    if any(strcmpi(fig_type,'fig'))
        saveas(i, [save_folder filesep save_name '_f_' num2str(i)], 'fig');
    end
    if any(strcmpi(fig_type,'png'))
%         saveas(i, [save_folder '\' save_name '_figure_' num2str(i)], 'png', '-r600');
        exportgraphics(figure(i),[save_folder filesep save_name '_figure_' num2str(i) '.png'],'Resolution',600)
    end
%   saveas(i, [save_folder '\' save_name '_figure_' num2str(i)], 'psc2')
%     if any(strcmpi(fig_type,'psc2'))
%         set(gcf, 'Renderer', 'painters');
%         saveas(i, [save_folder '\' save_name '_figure_' num2str(i)], 'psc2');
%     end
    if any(strcmpi(fig_type,'svg'))
        set(gcf, 'Renderer', 'painters');
        saveas(i, [save_folder filesep save_name '_figure_' num2str(i)], 'svg');
    end
end