%%% StreamFunction Vorticity
%%% LID DRIVEN CAVITY FLOW

close all
clear all
clc

%% Grid Parametres

np = 51;
dom = 1;
h = dom/(np-1);
x = 0:h:dom;
y = 0:h:dom;
Re= 100;
nu = 1/Re;
Ut = 1; % Lid Velocity

%% Initialization

psi = zeros(np,np);
w0 = zeros(np,np);

w_old = zeros(np,np);
u = zeros(np,np);      %Velocity in y direction
v = zeros(np,np);      %Velocity in x direction
u_sol = zeros(np,np);
v_sol = zeros(np,np);

% BC
w0(:,1) = -(2/(h*h))*(psi(:,2)); % left wall
w0(:,np) = -(2/(h*h))*(psi(:,np-1)); % right wall
w0(1,:) = -(2/(h*h))*(psi(2,:)) - (2*Ut/h); % top wall
w0(np,:) = -(2/(h*h))*(psi(np-1,:)); % bottom wall


%% Main loop

sol_index = 1;
t = 0;
dt = min(0.25*h*h/nu,4*nu/Ut*Ut); % Courant Number Formula
t_final = 30000*dt;

while t < t_final 
    
    % Boundary Conditions 
     
    w0(:,1) = -(2/(h*h))*(psi(:,2)); % left wall
    w0(:,np) = -(2/(h*h))*(psi(:,np-1)); % right wall
    w0(1,:) = -(2/(h*h))*(psi(2,:)) - (2*Ut/h); % top wall
    w0(np,:) = -(2/(h*h))*(psi(np-1,:)); % bottom wall

    % Vorticity Transport Equation
    
    w_old = w0;
    for i = 2:np-1
        for j = 2:np-1
    w0(i,j) = w_old(i,j) - (dt*((psi(i,j+1)-psi(i,j-1))/(2*h))*((w_old(i+1,j)-w_old(i-1,j))/(2*h)))...
        + (dt*((psi(i+1,j)-psi(i-1,j))/(2*h))*((w_old(i,j+1)-w_old(i,j-1))/(2*h)))...
        + (dt*nu*((w_old(i,j+1)+w_old(i,j-1)+w_old(i+1,j)+w_old(i-1,j) - 4*w_old(i,j))/(h*h)));
        end 
    end
    
    % Streamline Equation
    
    for i = 2:np-1
        for j =2:np-1
            psi(i,j)= ((h*h*w0(i,j))+psi(i,j+1)+psi(i,j-1)+psi(i+1,j)+psi(i-1,j))/4;   
        end
    end
 
    % Extracting Velocity
        
    v(1,1:np) =Ut; 
    for i = 2:np - 1
            for j = 2:np - 1
                u(i,j)=(psi(i,j+1)-psi(i,j-1))/(2*h);
                v(i,j)= -((psi(i+1,j)-psi(i-1,j))/(2*h));
            end
    end
    
    % Saving the solution after every dt
    for i = 1:np
        for j = 1:np
            u_sol(sol_index,i,j) = u(i,j);
            v_sol(sol_index,i,j) = v(i,j);
        end
    end
    sol_index = sol_index + 1;
    
    % Error Calculation
    
    error = 0;max_error = 1e-6;
    
    for i = 2:np-1
        for j =np-1
            error = error + abs(w0(i,j)-w_old(i,j))
        end
    end
    if error < max_error
        break;
    end
        
        t= t+dt;
end 

%% Plotting
    figure(1);
    x_dom = ((1:np)-1).*h;
    y_dom = 1-((1:np)-1).*h;
    [X,Y] = meshgrid(x_dom,y_dom);
    contourf(X,Y,v,40,'LineColor','none')
    colorbar
    colormap('hot')
    xlabel('x')
    ylabel('y')
    title(["StreamFunction Vorticity","Lid Driven Cavity flow"])
    figure(2);
    quiver(X,Y,v,u)

%% Animation

for i = 1:10:sol_index
    u_sample = v_sol(i,:,:);
    u_sample = reshape(u_sample, [size(u_sample,2) size(u_sample,3)]);
    x_dom = ((1:np)-1).*h;
    y_dom = 1-((1:np)-1).*h;
    [X,Y] = meshgrid(x_dom,y_dom);
    pcolor(X,Y,u_sample)
    shading interp
    colorbar
    colormap('jet')
    xlabel('x')
    ylabel('y')
    pause(0.25)
end