% For a rectangular grid

xgv = -2:0.1:2;
ygv = -2:0.1:2;
tv = 0:0.01:10;
dx = 0.1;
dy = 0.1;
dt = 0.01;

[X,Y] = meshgrid(xgv,ygv);

lx = length(xgv);
ly = length(ygv);
lt = length(tv);

phi = zeros(lx,ly,lt);
ux = zeros(lx,ly,lt);
uy = zeros(lx,ly,lt);

% Parameters
tau = 1;
epsilon = 0.1;
nu0 = 1;

% Initial conditions
for i=1:sx
    for j=1:sy
        if max(X(i,j),Y(i,j))>1 
            phi(i,j,1)=1;
        else phi(i,j,1)=-1;
        end
    end
end

% Edge velocities
ux(:,1,:) = u0;
uy(:,1,:) = 0;
ux(1,:,:) = 0;
uy(1,:,:) = -u0;
ux(:,ly,:) = -u0;
uy(:,ly,:) = 0;
ux(lx,:,:) = 0;
uy(lx,:,:) = u0;

% Corner velocities
ux(1,1,:) = 0;
uy(1,1,:) = 0;
ux(1,ly,:) = 0;
uy(1,ly,:) = 0;
ux(lx,1,:) = 0;
uy(lx,1,:) = 0;
ux(lx,ly,:) = 0;
uy(lx,ly,:) = 0;

ux_sym = sym(zeros(lx-2,ly-2));
for k = 1:numel(ux_sym)
    ux_sym(k) = sym(sprintf('ux%d', k));
end

uy_sym = sym(zeros(lx-2,ly-2));
for k = 1:numel(uy_sym)
    uy_sym(k) = sym(sprintf('uy%d', k));
end

str = ('initiator ');
for i=1:lx-2
    for j=1:ly-2
        eqn_x = ('');
        eqn_y = ('');
        str = strcat(str, eqn_x, eqn_y);
    end
end
str=strrep(str,'initiator ','');  



for k=1:lt
    for i=1:lx
        for j=1:ly
            if i==1 || i==lx || j==1 || j==ly
                phi(i,j,k+1) = 1;
            else
                phi(i,j,k+1) = phi(i,j,k) + dt*(-((ux(i,j,k)*((phi(i+dx,j,k)- ...
                    phi(i-dx,j,k))/(2*dx)))+(uy(i,j,k)*((phi(i,j+dy,k)-phi(i,j-dy,k))...
                    /(2*dy))))+(tau*((epsilon*(((phi(i+1,j,k)-2*phi(i,j,k)+ ...
                    phi(i-1,j,k))/((dx)^2))+((phi(i,j+1,k)-2*phi(i,j,k)+phi(i,j-1,k) ...
                    )/((dy)^2))))-((2*(phi(i,j,k))*(1-phi(i,j,k))*(1-2* ...
                    (phi(i,j,k))))/epsilon))));
            end
        end
    end
end




