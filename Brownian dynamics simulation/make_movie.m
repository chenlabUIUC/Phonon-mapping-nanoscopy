clear;

set(0, 'defaultTextInterpreter','latex');% set all text in figure as latex form
set(0, 'defaultLineLineWidth', 2);
set(0,'defaultTextFontSize',20);
set(0,'defaultAxesFontSize',20);
rng(1)
newcolors = rand([1000,1]);

 xmin = -100/17; xmax = 700/17;
 ymin = -100/17; ymax = 700/17;

v = VideoWriter('test_movie.avi');
v.FrameRate = 15;
open(v); 
M = moviein(10000);

FileNames = 'test_64particles';
k = 0;
% times = [0:100:3000,3400:400:5000,6000:1000:120000, 120000:200:150000];
% times = [0:100:3000,3400:400:5000,6000:1000:300000];
times = 1:1000:1300000;
    % figure/
for itime = 1:1:length(times)
    disp(itime)
    itime1 = times(itime);
    filename = [FileNames,'/dump_x.',num2str(itime1)];
        
    nskip = 9;
        
    fid=fopen(filename); 
    if fid == -1
        break;
    end
    datacell = textscan(fid, '%d %d %f %f %f', 'delimiter',' ', 'HeaderLines', nskip);
    fclose(fid);
    id_array = [datacell{2}];
    x_array = [datacell{3}];
    y_array = [datacell{4}]; %unit
    if itime == 1
        for i = 1:length(x_array)/4
            x1 = mean(x_array(4*i-3:4*i));
            y1 = mean(y_array(4*i-3:4*i));
            xx = 0.05*floor(x1/5);
            yy = 0.05*floor(y1/5);
            x_array(4*i-3:4*i) = x_array(4*i-3:4*i) - yy;
            y_array(4*i-3:4*i) = y_array(4*i-3:4*i) - xx;
        end
    end
    
    set(gcf,'unit','normalized','position',[0.1,0.2,0.3,0.6])
    
    % figure
    for i = 1:length(x_array)/4
        pgon  = polyshape(x_array(4*i-3:4*i),y_array(4*i-3:4*i));
        pg = plot(pgon);
        pg.FaceColor = newcolors(id_array(4*i):id_array(4*i)+3);
        hold on
    end
    hold off
    axis equal
    set(gca,'FontSize',26,'TickLabelInterpreter','latex')
    title(['timestep:',num2str(1000*itime-1000)]);
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    xticks([-300/17 -150/17 0 150/17 300/17, 450/17 600/17, 750/17 900/17])
    xticklabels({'-300','-150','0','150','300','450','600','750','900'})
    yticks([-300/17 -150/17 0 150/17 300/17, 450/17 600/17, 750/17 900/17])
    yticklabels({'-300','-150','0','150','300','450','600','750','900'})
    k = k + 1;
    M(:,k) = getframe(gcf);
    writeVideo(v,M(:,k));
end


% movie(M)
close(v)


X1 = datacell{1,3};
Y1 = datacell{1,4};
Xs = zeros(2,2);
for i = 1:2
    Xs(i,1) = mean(X1(i*4-3:i*4));
    Xs(i,2) = mean(Y1(i*4-3:i*4));
end

norm(Xs(1,:) - Xs(2,:))
