function [ result ] = center( image,x,y )
%CENTER star centroid
%   input:
%   image:normalization star field image
%   x:X-value of image coordinates
%   y:Y-value of image coordinates
%   output:
%   result:Center of mass of the star point
image=double(image);
result=zeros(length(x),2);
for i=1:length(result)
    sum_x=0;sum_y=0;s=0;
    for k=1:length(x{i})
        for h=1:length(y{i})
           s=s+image(x{i}(k),y{i}(h));
           sum_x=x{i}(k)*image(x{i}(k),y{i}(h))+sum_x;
           sum_y=y{i}(h)*image(x{i}(k),y{i}(h))+sum_y;
        end
    end
    result(i,1)=sum_x/s;
    result(i,2)=sum_y/s;
end
end


