function schematic(posx,posy)


%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
% idx = (rawStringColumns(:, 1) == "<undefined>");
% rawStringColumns(idx, 1) = "";
filename1 = 'PELOOP.00020000.dat';
%% Create output variable
fid = fopen(filename1);
dataall = textscan(fid,'%f %f %f %f %f %f','HeaderLines',1);
fclose(fid);
dataall = cell2mat(dataall);

xx = dataall(:,1); yy = dataall(:,2); zz = dataall(:,3);
uu = dataall(:,4); vv = dataall(:,5); ww = dataall(:,6);

%% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;
%% 对无法导入的数据进行的后处理。
% 在导入过程中未应用无法导入的数据的规则，因此不包括后处理代码。要生成适用于无法导入的数据的代码，请在文件中选择无法导入的元胞，然后重新生成脚本。

%% 创建输出变量

[nn,~]=size(xx);
% for jj=-10:4:10
%     z=x*0+jj;
%% for datacolor
xsize = max(xx(:)); ysize = max(yy(:)); zsize = max(zz(:));

xm = zeros(xsize,ysize,zsize); ym = xm; zm = xm;
um = xm; vm = xm; wm = xm;

eightcolor = [ ...
    200 100  17; ... 
    65  107 170; ...
    50  210   0; ...
    150 150 150; ...
    240 240 240; ...
    70  122   0; ...
    0     0 255; ...
    255 255   0; ...
    ]/255;
%     0 0 0; ... 
%     0 0 1; ...
%     0 1 0; ...
%     0 1 1; ...
%     1 0 0; ...
%     1 0 1; ...
%     1 1 0; ...
%     1 1 1; ...
for ii=1:nn
    xm(xx(ii),yy(ii),zz(ii)) = xx(ii);
    ym(xx(ii),yy(ii),zz(ii)) = yy(ii);
    zm(xx(ii),yy(ii),zz(ii)) = zz(ii);
    
%     deltaV = [uu(ii),vv(ii),ww(ii)];
        deltaV = [xx(ii)-32.5,yy(ii)-32.5,-zz(ii)];
    daltaVn = deltaV / (norm(deltaV)+eps);
    
    um(xx(ii),yy(ii),zz(ii)) = daltaVn(1);
    vm(xx(ii),yy(ii),zz(ii)) = daltaVn(2);
    wm(xx(ii),yy(ii),zz(ii)) = daltaVn(3);
end

boxones = ones(xsize,ysize,zsize);

nd = 1/sqrt(3);
datacolor111 = 1*um+1*vm+1*wm; datacolor111=nd*datacolor111;% 1  1  1
datacolor110 = 1*um+1*vm-1*wm; datacolor110=nd*datacolor110;% 1  1 -1
datacolor101 = 1*um-1*vm+1*wm; datacolor101=nd*datacolor101;% 1 -1  1
datacolor100 = 1*um-1*vm-1*wm; datacolor100=nd*datacolor100;% 1 -1 -1

datacolor011 =-1*um+1*vm+1*wm; datacolor011=nd*datacolor011;%-1  1  1
datacolor010 =-1*um+1*vm-1*wm; datacolor010=nd*datacolor010;%-1  1 -1
datacolor001 =-1*um-1*vm+1*wm; datacolor001=nd*datacolor001;%-1 -1  1
datacolor000 =-1*um-1*vm-1*wm; datacolor000=nd*datacolor000;%-1 -1 -1
%%
for ii=1:nn
    %     deltaV = [u(ii),v(ii),w(ii)]*0.1;
    %     datacolor = w(ii);
    %     datacolor = sqrt(x(ii)^2+y(ii)^2);
    if ( ...
            (xx(ii) == 1 || xx(ii)==64) ...
            ||(yy(ii)==1 || yy(ii)==64) ...
            ||(zz(ii)==16 ) ...
            ) ...
            && (mod(xx(ii),3)==1 && mod(yy(ii),3)==1 && mod(zz(ii),3)==1)
        %         if (xx(ii)==1&&xx(ii)==64&&mod(xx(ii)+1,3)==0 ...
        %             &&yy(ii)>=1&&yy(ii)<=64&&mod(yy(ii)+1,3)==0 ...
        %             &&(zz(ii)>=1&&zz(ii)<=15&&mod(zz(ii)+1,3)==0))
        
        disp([xx(ii),yy(ii),zz(ii)])
        %% 2019.04.19,15:30
        %         deltaV = [uu(ii),vv(ii),ww(ii)]*0.8;
        %         pos=[xx(ii),yy(ii),zz(ii)]-0.0*deltaV;
        %         datacolor=norm(deltaV);
        % %%
        %         arrow3D(pos,deltaV,ww(ii));
        %%
%         deltaV = [uu(ii),vv(ii),ww(ii)]*1.5;
                deltaV = [xx(ii)-32.5,yy(ii)-32.5,-zz(ii)];
        deltaV = 8*deltaV / (norm(deltaV)+eps);
        
        pos=[xx(ii)+posx,yy(ii)+posy,zz(ii)]-0.5*deltaV;
        
        datacolorp = [datacolor111(xx(ii),yy(ii),zz(ii));...
            datacolor110(xx(ii),yy(ii),zz(ii));...
            datacolor101(xx(ii),yy(ii),zz(ii));...
            datacolor100(xx(ii),yy(ii),zz(ii));...
            datacolor011(xx(ii),yy(ii),zz(ii));...
            datacolor010(xx(ii),yy(ii),zz(ii));...
            datacolor001(xx(ii),yy(ii),zz(ii));...
            datacolor000(xx(ii),yy(ii),zz(ii))];
        
        
        datacolor = find(datacolorp==max(datacolorp),1);
        %%
        if norm(deltaV)>0
            arrow3D(pos,deltaV,datacolor);
        end
    end
    
end
%%
hold on
x=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0];
y=[0 0 1 0 0 0;0 1 1 1 0 0;0 1 1 1 1 1;0 0 1 0 1 1];
z=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1];

xxb = x*68+posx;
yyb = y*68+posy;
zzb = z*14; zzb = zzb+0;
hbox1=fill3(xxb,yyb,zzb, 'r','facealpha',0.4,'FaceColor',[0,0,0],'EdgeAlpha','0.8');
%%
hold on

xxb = x*60+2+posx;
yyb = y*60+2+posy;
zzb = z*7; zzb = zzb+7.5;
% hbox2=fill3(xxb,yyb,zzb, 'r','facealpha',0.7,'FaceColor',[0,0,0],'EdgeAlpha','0.7');
%%
axis equal tight off
colormap(eightcolor)
view(-102,36)

caxis([0.5,8.5])

h=light;
lightangle(h,45,45)
h1=light;
lightangle(h1,-45,-45)
lighting gouraud;
