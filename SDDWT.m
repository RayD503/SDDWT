function [starCluster] = SDDWT(img_,T)
%SDDWT Star detection based on dyadic wavelet transform
%   Detailed explanation goes here
%   input:
%   img_:normalization star field image
%   T:threshold
%   output
%   starCluster:Collection of star points

lodyadf0=sqrt(2).*[0.125 0.375 0.375 0.125];
hidyadf0=sqrt(2).*[-0.25 0.5 -0.25];
L=2; %Decomposition scale

tic
dwt = zeros(size(img_));
for i=1:length(img_)
    lodyadf=lodyadf0;
    hidyadf=hidyadf0;
    dwt(:,i) = img_(i,:)';
    for d = 1:L-1 %Low frequency section
        s = dwt(:,i)';
        dwt(:,i) = iconv(lodyadf,s)';
        
        f = zeros(1,2*length(lodyadf));
        f(1:2:2*length(lodyadf)-1) = lodyadf;%Low frequency filter update
        
        f2 = zeros(1,2*length(hidyadf));
        f2(1:2:2*length(hidyadf)-1) = hidyadf;%high frequency filter update
        
        lodyadf = f;
        hidyadf = f2;
    end
    s = dwt(:,i)';
    dwt(:,i) = iconv(hidyadf,s)';%high frequency section
    for j = 1:floor(1.25*(2^(L)-1))
        p = lshift(dwt(:,i)');
        dwt(:,i) = p';
    end
end

% maxima
dwt=dwt';
[rows,cols]=size(dwt);
localmaxima_row = zeros(rows,cols);
localmaxima_row_all= zeros(rows,cols);
t      = 1:rows;
tplus  = rshift(t);
tminus = lshift(t);

for k=1:rows
    x=(dwt(k,:));
    x = max([x(t); x(tplus); x(tminus)]);
    localmaxima_row(k,:) = (dwt(k,:))>=x;
    localmaxima_row_all(k,:)=localmaxima_row(k,:).*dwt(k,:);
    localmaxima_row(k,:) = localmaxima_row(k,:) .* ((dwt(k,:))>T);
end
x=[];
toc
localmaxima_row_=localmaxima_row;
%% Separate individual star points
x=[];y=[];
n=1;%Number of stars
Ltemp=3;%L
v=1;%CARµÄÆ«ÒÆ
for i=4:length(localmaxima_row(:,1))-3
    for k=4:length(localmaxima_row(1,:))-3
        tl=1;
        if localmaxima_row(i,k)==1%Get CAR
            x1=i;%Upper boundary
            x2=x1;%Lower boundary
            localmaxima_row(i,k)=0;
            k_=k;
            
            while 1
                if x1+tl>length(localmaxima_row(:,1))-3 
                    break;
                end    
                [~,lib]=find(localmaxima_row(x1+tl,k_(end)-v:k_(end)+v)==1);
                if isempty(lib)
                    if tl==1
                        break;
                    end
                    
                    kmin=min(k_);kmax=max(k_);
                    if x1-Ltemp>0 && x2+Ltemp<1025 && kmin-round(tl/2)-Ltemp>0 && kmax+round(tl/2)+Ltemp<1025
                        starCluster(n,:)=[x1-Ltemp,x2+Ltemp,kmin-round(tl/2)-Ltemp,kmax+round(tl/2)+Ltemp-1];
                        localmaxima_row(x1-Ltemp:x2+Ltemp,kmin-round(tl/2)-Ltemp:kmax+round(tl/2)+Ltemp-1)=0;
                        n=n+1;
                        break;
                    else
                        break;
                    end
                else
                    x2=x1+tl;
                    tl=tl+1;
                    k_=[k_;k_(end)-v-1+lib(1)];
                end
            end
        else
            continue;
        end
        k_=[];
    end
end
end

