function Plot_Rotary(Fig,Ax,NP,r,alpha,phi,phi_s,phi_ss,phi_sss,varphi,varphi_s,varphi_ss,V_alpha,V_alpha_s,V_alpha_ss)

L_Frame = 0.5/NP(1);
% Compute Frenet frame at the nodes
[r_NC r_NC_s r_NC_ss r_NC_sss alpha_NC alpha_NC_s alpha_NC_ss ] = Compute_Functions_NC(NP,r,alpha,phi,phi_s,phi_ss,phi_sss,varphi,varphi_s,varphi_ss,V_alpha,V_alpha_s,V_alpha_ss); 

% Find the values just at the nodes
r_s  = r_NC_s(1:NP(2):end,:);
r_ss = r_NC_ss(1:NP(2):end,:);

Vec_t  = zeros(NP(1)+1,3);
Vec_n  = zeros(NP(1)+1,3);
Vec_b  = zeros(NP(1)+1,3);
Vec_d1 = zeros(NP(1)+1,3);
Vec_d2 = zeros(NP(1)+1,3);
kappa  = zeros(NP(1)+1,1);
for i=1:1:NP(1)+1
    % Tangent vector
    Vec_t(i,:) = r_s(i,:)/sqrt(r_s(i,:)*r_s(i,:)');
    % Curvature
    kappa(i)   = sqrt(r_ss(i,:)*r_ss(i,:)');
    % Normal vector
    if abs(kappa(i))>0
       % AUX        = cross(r_ss(i,:),r_s(i,:));
       % Vec_n(i,:) = cross(r_s(i,:),AUX);
       % Vec_n(i,:) =  Vec_n(i,:)/sqrt( Vec_n(i,:)* Vec_n(i,:)');  
       Vec_n(i,:) =  r_ss(i,:)/kappa(i);
    end
    % Binormal vector
    Vec_b(i,:) = cross(squeeze(Vec_t(i,:)),squeeze(Vec_n(i,:)));
    % Vector d1
    Vec_d1(i,:)    = cos(alpha(i))*Vec_n(i,:) - sin(alpha(i))*Vec_b(i,:);
    % Vector d2
    Vec_d2(i,:)    = sin(alpha(i))*Vec_n(i,:) + cos(alpha(i))*Vec_b(i,:);
end

set(Fig,'CurrentAxes',Ax)
cla(Ax)
plot3(r(:,1),r(:,2),r(:,3),'.r')
hold on
plot3(r(:,1),r(:,2),r(:,3),'-r')

% Plot Frenet frame
for i=1:1:NP(1)+1
     plot3([r(i,1) r(i,1)+L_Frame*Vec_t(i,1)] ,[r(i,2) r(i,2)+L_Frame*Vec_t(i,2)] ,[r(i,3) r(i,3)+L_Frame*Vec_t(i,3)],'b')
     plot3([r(i,1) r(i,1)+L_Frame*Vec_d1(i,1)],[r(i,2) r(i,2)+L_Frame*Vec_d1(i,2)],[r(i,3) r(i,3)+L_Frame*Vec_d1(i,3)],'b')
     plot3([r(i,1) r(i,1)+L_Frame*Vec_d2(i,1)],[r(i,2) r(i,2)+L_Frame*Vec_d2(i,2)],[r(i,3) r(i,3)+L_Frame*Vec_d2(i,3)],'b')
end

set(Ax,'ZDir','Reverse');
set(Ax,'XDir','Reverse');
xlabel('x')
ylabel('y')
zlabel('z')
view(-37.5, 30)
grid on
axis([-1 0 -0.5 0.5 -1 0])
pause(0.1)

end