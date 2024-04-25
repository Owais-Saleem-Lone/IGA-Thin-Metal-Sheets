% Get the location of control points where force or DCs are applied

[xs,ys,zs] = createSupportsForIGABernoulliBeam2D(CP,homDOFs);
[xf,yf,zf] = createForceArrowsForIGABernoulliBeam2D(CP,F);

%% 5. Plot the supports on the geometry
for k =1:length(xs(:,1))
    plot3(app.UIAxes2_3,xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
%% 6. Plot the load arrows on the geometry
for k =1:length(xf(:,1))
    plot3(app.UIAxes2_3,xf(k,:),yf(k,:),zf(k,:),'Color','blue','Linewidth',5);
    plot3(app.UIAxes2_3,xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',5);
end
%%
function [xs, ys, zs] = createSupportsForIGABernoulliBeam2D ...
    (CP, rb)
nu = length(CP(:,1));
% scaling factors for the support triangles
maximum = max(max(max(max(CP))));
minimum = min(min(min(min(CP))));

% Average the factor with respect to the maximum and minimum values
factor = (maximum-minimum)/5;

% Initialize the output arrays
xs = zeros(length(rb),4);
ys = zeros(length(rb),4);
zs = zeros(length(rb),4);

for k = 1:length(rb)
    % Get the corresponding Control Point number p and indices CP(i,j)
    h=rb(k)/2;
    p=ceil(h);
    j=ceil(p/nu);
    i=p-(j-1)*nu;

    %(rb is odd -> horizontal support)
    if (p~=h)
        xs(k,1)=CP(i,1);
        xs(k,2)=CP(i,1)-0.1732*factor;
        xs(k,3)=CP(i,1)-0.1732*factor;
        xs(k,4)=xs(k,1);
        ys(k,1)=CP(i,2);
        ys(k,2)=CP(i,2)+0.1*factor;
        ys(k,3)=CP(i,2)-0.1*factor;
        ys(k,4)=ys(k,1);
        zs(k,1:4)=CP(i,3);
        %(rb is even -> vertical support)
    else
        xs(k,1)=CP(i,1);
        xs(k,2)=CP(i,1)-0.1*factor;
        xs(k,3)=CP(i,1)+0.1*factor;
        xs(k,4)=xs(k,1);
        ys(k,1)=CP(i,2);
        ys(k,2)=CP(i,2)-0.1732*factor;
        ys(k,3)=CP(i,2)-0.1732*factor;
        ys(k,4)=ys(k,1);
        zs(k,1:4)=CP(i,3);
    end
end

end
%%

function [xf, yf, zf] = createForceArrowsForIGABernoulliBeam2D ...
    (CP, F)
nf = sum(F~=0);
xf = zeros(nf,2);
yf = zeros(nf,2);
zf = zeros(nf,2);
% Initialize counters
k = 1;
l = 1;
% Loop over all the Control Points
for i = 1:length(CP(:,1,1))
    if F(k)~=0
        xf(l,1) = CP(i,1);
        xf(l,2) = CP(i,1)-F(k)/max(abs(F));
        yf(l,1) = CP(i,2);
        yf(l,2) = CP(i,2);
        zf(l,1) = CP(i,3);
        zf(l,2) = CP(i,3);

        % Update internal counter
        l = l + 1;
    end
    if F(k+1)~=0
        xf(l,1) = CP(i,1);
        xf(l,2) = CP(i,1);
        yf(l,1) = CP(i,2);
        yf(l,2) = CP(i,2)-F(k+1)/max(abs(F));
        zf(l,1) = CP(i,3);
        zf(l,2) = CP(i,3);
        l = l + 1;
    end
    k = k + 2;
end
end