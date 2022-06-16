clc
clear
close all

%% Generate uneven backgrounds
img=ones(1024,1024)*30;
sigma=1000;
X=1:1024;
Y=1:1024;
[XX,YY]=meshgrid(X,Y);
Z=(XX-512).^2+(YY-512).^2;
Z=-Z/(2*sigma^2);
Z=exp(Z)/(2*pi*sigma^2);

G_min=min(min(Z));
G_max=max(max(Z));
Z=(Z-G_min)./(G_max-G_min)*30;

img=Z+img;
%% Randomly distributed stars
load 'RandomlyDistributedPoints.mat'

Mv=1:0.05:5.95;
simga=1.5;
A=255-30*Mv;

for i=1:100
    star{i}=(XX-coor(i,1)).^2+(YY-coor(i,2)).^2;
    star{i}=-star{i}./(2*simga^2);
    star{i}=A(i).*exp(star{i});
end

for i=1:100
    img=max(img,star{i});
end
%% white Gauss noise
dB=10;%Power of nosie
w=wgn(1024,1024,dB);

img=img+w;

T=10^(dB/20)/133.5;% threshold T
disp('The threshold:')
disp(T)
%% Show the simulated star field image
figure
surf(img,'FaceColor','interp','EdgeColor','none')
xlabel('X axis')
ylabel('Y axis')
zlabel('Gray value')
colorbar
set(gca,'Fontsize',20)
axis([1 1024 1 1024 0 255])

c=colorbar;
ax=gca;
axpos=ax.Position;
c.Position(3)=0.5*c.Position(3);
ax.Position=axpos;

%% SDDWT
%normalization
img_=img/255;

starCluster = SDDWT(img_,T);
%% Star centroid
for i=1:length(starCluster(:,1))
    x{i}=starCluster(i,1):starCluster(i,2);
    y{i}=starCluster(i,3):starCluster(i,4);
end

point= center( img_,x,y );
%% Detection rate
k=1;
for i=1:length(point)
    temp=sqrt((point(i,1)-coor(:,2)).^2+(point(i,2)-coor(:,1)).^2);
    [a,b]=min(temp);
    if a<3
        point_(k,1:2)=point(i,1:2);
        point_(k,3)=b;
        k=k+1;
    end
end

disp('real stars detected')
disp(length(point_))
disp('false stars detected')
disp(length(point)-length(point_))

