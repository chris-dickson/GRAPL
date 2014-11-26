function grid_Movie(fn_prefix,no_frames,outfile,width,height,fps)

    % color map lookup function
    function [r,g,b] = color_lookup(x)
        if (x <= 1)
            r = 1-x; g = 1-x; b = 1;
        else
            if (x > 4)
                r = 1; g = 0; b = 0;
            else
                r = 1; g = 1-x/4; b = 0;
            end
        end
    end


    if (no_frames == 0)
        
        %read in output
        filename = strcat(fn_prefix, strcat('_0.txt'));
        fp = fopen(filename,'r');
        while feof(fp) == 0
            x = fscanf(fp,'%d',1);
            x = x+1;

            y = fscanf(fp,'%d',1);
            y = y+1;

            d = fscanf(fp,'%g',1);

            if ( ~isempty(x) && ~isempty(y) && ~isempty(d) )
                Z(y,x) = d;
            end
        end
        status = fclose(fp); 
        
        
        
        for i = 1:2*height-1
           for j = 1:2*width-1
               [r,g,b] = color_lookup(Z(i,j));
               
               C(i,j,1) = r;
               C(i,j,2) = g;
               C(i,j,3) = b;
           end
       end
        
       image(C);
       axis([0 2*width 0 2*height]);
       set(gcf,'Position',[10 10 1000 600]);
       
       return;
    end

    aviobj = avifile(outfile,'fps',fps);

 
    for l = 0:no_frames

        filename = strcat(fn_prefix, strcat('_', strcat(int2str(l),'.txt')))

        [X,Y] = meshgrid(0:2*(width-1),0:2*(height-1));

        X = X + ones(2*height-1,2*width-1);
        Y = Y + ones(2*height-1,2*width-1);

        Z = zeros(2*height-1,2*width-1);

        %read in output
        fp = fopen(filename,'r');
        while feof(fp) == 0
            x = fscanf(fp,'%d',1);
            x = x+1;

            y = fscanf(fp,'%d',1);
            y = y+1;

            d = fscanf(fp,'%g',1);

            if ( ~isempty(x) && ~isempty(y) && ~isempty(d) )
                Z(y,x) = d;
            end
        end
        status = fclose(fp);
        
       

        %title_ext = strcat(' :  ', strcat(int2str(l),strcat('/',int2str(no_frames))));        
        %grid_title = strcat(title_prefix,title_ext);
        
       for i = 1:2*height-1
           for j = 1:2*width-1
               [r,g,b] = color_lookup(Z(i,j));
               
               C(i,j,1) = r;
               C(i,j,2) = g;
               C(i,j,3) = b;
           end
       end
        
        
%        contour3(X,Y,Z);
%        surface(X,Y,Z);
         image(C);

     
   %     colormapeditor;       

    %    mycmap = get(gcf,'Colormap');
    %    save('MyColormaps','mycmap');
       
        
        axis([0 2*width 0 2*height]);
   %     axis manual;
%        title(grid_title);
        set(gcf,'Position',[10 10 1000 600]);
        
        
        one_frame = getframe(gca);    
        aviobj = addframe(aviobj,one_frame);
              
    end
    aviobj = close(aviobj);   
end