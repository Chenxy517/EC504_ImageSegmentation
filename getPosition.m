
function getPosition()
    close all; clear all;clc;
    % 创建figure和按钮
    f = figure('Units','normalized','Position',[0.2 0.2 0.5 0.5]);
    uicontrol('Style','pushbutton','String','Select Image','Position',[20 20 100 30],'Callback',@select_image);
    
    % 创建两个text box
    uicontrol('Style','text','String','X:','Position',[150 20 20 30],'HorizontalAlignment','right');
    edit1 = uicontrol('Style','edit','Position',[175 25 50 20],'HorizontalAlignment','left');
    uicontrol('Style','text','String','Y:','Position',[250 20 20 30],'HorizontalAlignment','right');
    edit2 = uicontrol('Style','edit','Position',[275 25 50 20],'HorizontalAlignment','left');
    
    % 显示图片的axes
    axes_handle = axes('Parent',f,'Units','normalized','Position',[0.2 0.3 0.6 0.6]);
    axis(axes_handle, "square");

    % 选择图片的回调函数
    function select_image(~,~)
        [filename,pathname] = uigetfile({'*.jpg;*.png;*.bmp;*.gif;*.tif','Image Files (*.jpg,*.png,*.bmp,*.gif,*.tif)';'*.*','All Files (*.*)'});
        if isequal(filename,0) || isequal(pathname,0)
            return;
        end
        
        img = imread(fullfile(pathname,filename));
        imshow(imresize(img, [600, 600]), 'Parent', axes_handle);

        % 新的axes图层
        ax = axes('Parent', f, 'Position', get(axes_handle,'Position'));
        axis(ax, "square");
        set(ax, 'Color', 'none');
        set(ax,'ButtonDownFcn',@record_coords);
        
    end



    % 记录坐标的回调函数
    function record_coords(src,event)
        % 获取鼠标点击的坐标
        x = 1;
        point = get(src,'CurrentPoint');
        x = point(1,1);
        y = point(1,2);
        fprintf('Node x=%f,y=%f\n',x,y);
        
        % 显示坐标
        set(edit1,'String',num2str(x));
        set(edit2,'String',num2str(y));

    end
end


