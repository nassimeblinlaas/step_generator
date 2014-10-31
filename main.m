%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hirukawa multicontact step generator 2007
%
% Implémenté par Nassime BLIN et Maximilien Naveau 2014
% pour Gepetto LAAS - CNRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all ;
clc ;
clear ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delarations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global uLINK G
G = 9.81 ;

number_of_samples = 100; % size of data to treat
pZ = 0.6487;            % position Z of robot constant
period = 0.005;         % sampling period in seconds

WAIST = 1;              % labels
RLEG_JOINT0 = 2;        %
RLEG_JOINT1 = 3;        %
RLEG_JOINT2 = 4;        %
RLEG_JOINT3 = 5;        %
RLEG_JOINT4 = 6;        %
RLEG_JOINT5 = 7;        %
LLEG_JOINT0 = 8;        %
LLEG_JOINT1 = 9;        %
LLEG_JOINT2 = 10;       %
LLEG_JOINT3 = 11;       %
LLEG_JOINT4 = 12;       %
LLEG_JOINT5 = 13;       %
CHEST_JOINT0 = 14;      %
CHEST_JOINT1 = 15;      %
HEAD_JOINT0 = 16;       %
HEAD_JOINT1 = 17;       %
RARM_JOINT0 = 18;       %
RARM_JOINT1 = 19;       %
RARM_JOINT2 = 20;       %
RARM_JOINT3 = 21;       %
RARM_JOINT4 = 22;       %
RARM_JOINT5 = 23;       %
RARM_JOINT6 = 24;       %
LARM_JOINT0 = 25;       %
LARM_JOINT1 = 26;       %
LARM_JOINT2 = 27;       %
LARM_JOINT3 = 28;       %
LARM_JOINT4 = 29;       %
LARM_JOINT5 = 30;       %
LARM_JOINT6 = 31;       %

v_ref_B  = zeros(3,1);  % waist linear speed
v_ref_F_1 = zeros(3,1);
v_ref_F_2 = zeros(3,1);
v_ref_H_1 = zeros(3,1);
v_ref_H_2 = zeros(3,1);

w_ref_B  = zeros(3,1);
w_ref_F_1 = zeros(3,1);
w_ref_F_2 = zeros(3,1);
w_ref_H_1 = zeros(3,1);
w_ref_H_2 = zeros(3,1);

xi_B  = zeros(6,1);
xi_F_1 = zeros(6,1);
xi_F_2 = zeros(6,1);
xi_H_1 = zeros(6,1);
xi_H_2 = zeros(6,1);

d_theta_leg_1 = zeros(6,1);
d_theta_leg_2 = zeros(6,1);
d_theta_arm_1 = zeros(6,1);
d_theta_arm_2 = zeros(6,1);

q_values = zeros(number_of_samples, 30);

I3 = eye(3,3);

sample = 0;             % current sample

L_prev = zeros(3,1);    % L value of n-1 sample
x_G_prev = 0;
y_G_prev = 0;
err = [0.1; 0.1];       % converging threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Initialisation\n')

uLINK = loadHRPdata('HRP2main_full.wrl');

fprintf('Reading ./TestMorisawa2007ShortWalk32TestFGPI.dat\n')
%Whole_data = csvread('./morisawa.csv');
Whole_data = load('./TestMorisawa2007ShortWalk32TestFGPI.dat');
Data = Whole_data(1:number_of_samples, :);
given_xi_B = zeros(size(Data), 6);
res_xi_B = zeros(size(Data), 6);
halfsitting = load('./halfsitting.dat');
fprintf('Reading done\n')

q_values(1, :) = halfsitting * pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build contact points data
is_contact = zeros(length(Data) , 2);
contact_coord = zeros(length(Data), 4);
% three scalars per vector foreach two contact points
normal_vectors = zeros(length(Data), 3 * 2);    

for i = 1:number_of_samples
    if Data(i, 13) == 0
       is_contact(i, 1) = 1;
       contact_coord(i, 1) = Data(i, 11);       % coord x LF
       contact_coord(i, 2) = Data(i, 12);       % coord y LF
       contact_coord(i, 3) = Data(i, 13);       % coord z LF
       normal_vectors(i, 1) = 0;                % x1
       normal_vectors(i, 2) = 0;                % y1
       normal_vectors(i, 3) = 1;                % z1
    end
    if Data(i, 25) == 0
        is_contact(i, 2) = 1;
        contact_coord(i, 4) = Data(i, 23);      % coord x RF
        contact_coord(i, 5) = Data(i, 24);      % coord y RF
        contact_coord(i, 6) = Data(i, 25);      % coord z RF
        normal_vectors(i, 4) = 0;               % x2
        normal_vectors(i, 5) = 0;               % y2
        normal_vectors(i, 6) = 1;               % z2
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = TotalMass(1);                               % total robot mass

v_CoM = [ Data(sample, 6); Data(sample, 7); 0];

uLINK(WAIST).p = [ 0 0 0.6487 ]' ;
uLINK(WAIST).R = eye(3,3) ;
uLINK(WAIST).v = v_CoM + cross(r_bc,w_ref_B) ;
uLINK(WAIST).w = [ 0; 0; 0] ;
for i = 2:length(uLINK)
    uLINK(i).q   = halfsitting(i-1) * pi/180;
    uLINK(i).dq  = 0.0 ;
    uLINK(i).ddq = 0.0 ;
end

ForwardKinematics(1);

for i = 1:length(uLINK)
    uLINK(i).v  = [0.0 ; 0.0 ; 0.0 ] ;
    uLINK(i).w  = [0.0 ; 0.0 ; 0.0 ] ;
end

ForwardVelocity(1);

CoM_init = calcCoM();
r_bc = uLINK(WAIST).p - CoM_init;     % vector from CoM to Base Link origin


 hold off
    newplot
    DrawAllJoints(1);
    axis equal
  	set(gca,...
        'CameraPositionMode','manual',...
        'CameraPosition',[4,4,1],...
        'CameraViewAngleMode','manual',....
        'CameraViewAngle',15,...
        'Projection','perspective',... 
        'XLimMode','manual',...
        'XLim',[-0.5 0.5],...
        'YLimMode','manual',...
        'YLim',[-0.5 0.5],...
        'ZLimMode','manual',...
        'ZLim',[0 1.5])
    grid on
    text(0.5, -0.4, 1.4, ['time=',num2str(sample*period,'%5.3f')])
    drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Big loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sample = 1 : number_of_samples
    if (mod(sample, 10) == 0) sample   % for debug usage
    end
    
    ForwardKinematics(1);
    ForwardVelocity(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1 : give waist linear and angular speed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v_CoM = [ Data(sample, 6); Data(sample, 7); 0];
    
    v_ref_B   = v_CoM + cross(r_bc,w_ref_B) ;     % waist speed vector
    w_ref_B   = [ 0; 0; 0];
    
    uLINK(WAIST).v = v_ref_B ;
    
    v_ref_F_1 = [ Data(sample, 14); Data(sample, 15); Data(sample, 16)];
    w_ref_F_1 = [ 0; 0; 0];
    
    v_ref_F_2 = [ Data(sample, 26); Data(sample, 27); Data(sample, 28)];
    w_ref_F_2 = [ 0; 0; 0];
    
    v_ref_H_1 = [ 0; 0; 0];
    w_ref_H_1 = [ 0; 0; 0];
    
    v_ref_H_2 = [ 0; 0; 0];
    w_ref_H_2 = [ 0; 0; 0];
    
    %[0.1 ; 0.1 ; 0.1];


    iteration   = 0;                            % algorithm iteration
    converge    = 0;                            % convergence boolean


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Small loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %while ( converge == 0 )
    while ( iteration < 3 )
        iteration = iteration + 1;

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 2 : find bodies angular speeds 
        % having new linear and angular speed of B
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % vector waist -> end effector
        r_B_F1 = hat(uLINK(WAIST).p - uLINK(RLEG_JOINT5).p);
        r_B_F2 = hat(uLINK(WAIST).p - uLINK(LLEG_JOINT5).p);
        r_B_H1 = hat(uLINK(WAIST).p - uLINK(RARM_JOINT5).p);
        r_B_H2 = hat(uLINK(WAIST).p - uLINK(LARM_JOINT5).p);

        xi_B   = [ v_ref_B   ; w_ref_B   ];
        xi_F_1 = [ v_ref_F_1 ; w_ref_F_1 ];
        xi_F_2 = [ v_ref_F_2 ; w_ref_F_2 ];
        xi_H_1 = [ v_ref_H_1 ; w_ref_H_1 ];
        xi_H_2 = [ v_ref_H_2 ; w_ref_H_2 ];

        % leg 1, find angular speeds d_theta
        tmp = [eye(3,3),-r_B_F1;zeros(3,3),eye(3,3)];
        route_leg1 = FindRoute(RLEG_JOINT5);
        J_leg_1 = CalcJacobian(route_leg1);
        d_theta_leg_1 = J_leg_1\ xi_F_1 - J_leg_1\ tmp * xi_B;

        % leg 2
        tmp = [eye(3,3),-r_B_F2;zeros(3,3),eye(3,3)];
        route_leg2 = FindRoute(LLEG_JOINT5);
        J_leg_2 = CalcJacobian(route_leg2);
        d_theta_leg_2 = J_leg_2\ xi_F_2 - J_leg_2\ tmp * xi_B;

        % arm 1
        tmp = [eye(3,3),-r_B_H1;zeros(3,3),eye(3,3)];
        route_arm1 = FindRoute(RARM_JOINT5);
        route_arm1 = route_arm1(3:end);
        J_arm_1 = CalcJacobian(route_arm1);
        d_theta_arm_1 = J_arm_1\ xi_H_1 - J_arm_1\ tmp * xi_B;

        % arm 2
        tmp = [eye(3,3),-r_B_H2;zeros(3,3),eye(3,3)];
        route_arm2 = FindRoute(LARM_JOINT5);
        route_arm2 = route_arm2(3:end);
        J_arm_2 = CalcJacobian(route_arm2);
        d_theta_arm_2 = J_arm_2\ xi_H_2 - J_arm_2\ tmp * xi_B;

        for i = 1:length(uLINK)
            uLINK(i).dq = 0.0 ;
        end
        for i = 1:length(route_arm2)
            uLINK(route_leg1(i)).dq = d_theta_leg_1(i);
            uLINK(route_leg2(i)).dq = d_theta_leg_2(i);
            uLINK(route_arm1(i)).dq = d_theta_arm_1(i);
            uLINK(route_arm2(i)).dq = d_theta_arm_2(i);
        end
        ForwardKinematics(1);
        ForwardVelocity(1);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3 : find angular momentum L
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Call MH subroutine finding intertia matrices
        [ M_d_theta, H_d_theta, I_tilde, M_leg_1, M_leg_2, ...
            M_arm_1, M_arm_2, H_leg_1, H_leg_2, H_arm_1, H_arm_2 ] = MH();
        d_theta = [d_theta_leg_1' d_theta_leg_2' ...
            d_theta_arm_1' d_theta_arm_2'];

        A = [ M*I3 , -M*hat(r_bc) , M_d_theta ;
              zeros(3,3) , I_tilde , H_d_theta ];

        B = [ v_ref_B ; w_ref_B ; d_theta' ]  ;

        PL = A * B;
        P = PL(1:3);                            % linear momentum
        L = PL(4:6);                            % angular momentum

        P_kajita(sample,:) = calcP(1);
        L_kajita(sample,:) = calcL(1);
    
        P_nous(sample,:) = P;
        L_nous(sample,:) = L;
        
        dL = (L - L_prev) / period;             % finite difference method

        L_prev = L;                             % save previous value

        [M2,mc2,c2,I_tilde2] = calcMHzero(1);
        
        
        M2_leg_1 = [uLINK(route_leg1).mj] ; H2_0_leg_1 = [uLINK(route_leg1).hj]
        M2_leg_2 = [uLINK(route_leg2).mj] ; H2_0_leg_2 = [uLINK(route_leg2).hj]
        M2_arm_1 = [uLINK(route_arm1).mj] ; H2_0_arm_1 = [uLINK(route_arm1).hj]
        M2_arm_2 = [uLINK(route_arm2).mj] ; H2_0_arm_2 = [uLINK(route_arm2).hj]
        
        M2_d_theta   = [M2_leg_1  , M2_leg_2  , M2_arm_1  , M2_arm_2  ];
        H2_0_d_theta = [H2_0_leg_1, H2_0_leg_2, H2_0_arm_1, H2_0_arm_2];
        H2_d_theta   = H2_0_d_theta - hat(uLINK(WAIST).c_tilde) * M2_d_theta ; 
        
        A2 = [ M*I3 , -M*hat(r_bc) , M2_d_theta ;
                  zeros(3,3) , I_tilde2 , H2_d_theta ];

        B2 = [ v_ref_B ; w_ref_B ; d_theta' ]  ;

        PL2 = A * B;
        P2 = PL2(1:3);                            % linear momentum
        L2 = PL2(4:6);                            % angular momentum

        P
        P2
        calcP = calcP(1)
        
        L
        L2
        calcL = calcL(1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 4 : give \ddot{z_G}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % suppose that :
        ddZg = 0; % always



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 5 : give lambda_k, slope of the ground
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % hypothesis : slope is zero
        % => alpha = 0   -> equations (16) (17) (18)
        alpha = 0;

        % for two contact points :
        lambda = zeros(2, 1);
        lambda(1) = 0.5;
        lambda(2) = 0.5;

        epsilon = zeros(2, 1);
        epsilon_ref = zeros(2, 1);

        % the multiplication by 1 is the normal vector n_k_z
        epsilon(1) = (1 - alpha) * M * ( ddZg + G ) * ...
            ( ( lambda(1) ) / ( lambda(1) * 1 + lambda(2) * 1 ) ); 
        epsilon(2) = (1 - alpha) * M * ( ddZg + G ) * ...
            ( ( lambda(2) ) / ( lambda(2) * 1 + lambda(2) * 1 ) ); 


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 6 : find d_d_x_G and d_d_y_G
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Init
        tau_C_x = 0;
        tau_C_y = 0;
        K = length(epsilon);
        C =   zeros(length(Data), 3);           % contact point x, y, z 
        x_G = zeros(length(Data), 3);           % for x_G \dotx_G \ddotx_G
        y_G = zeros(length(Data), 3);           % for y_G \doty_G \ddoty_G


        % Torques
        for k = 1:K
            offset = (k-1) * 3;
            x = offset + 1;
            y = offset + 2;
            z = offset + 3;

            if ( is_contact(sample, k) == 1 )
                tau_C_x = tau_C_x + epsilon(k) * ...
                    ( contact_coord(sample, y) * ...
                    normal_vectors(sample, z) - ...
                    contact_coord(sample, z) * normal_vectors(sample, y) );

                tau_C_y = tau_C_y + epsilon(k) * ...
                    ( contact_coord(sample, x) * ...
                    normal_vectors(sample, z) - ...
                    contact_coord(sample, z) * normal_vectors(sample, x) );
            end
        end 


        % Contact points
        epsilon_tot = sum(epsilon);
        for k = 1:K
            offset = (k-1) * 3;
            x = offset + 1;
            y = offset + 2;
            z = offset + 3;

            C(sample, 1) = C(sample, 1) + ( epsilon(k) / epsilon_tot) * ...
                contact_coord(sample, x);   % x_C
            C(sample, 2) = C(sample, 2) + ( epsilon(k) / epsilon_tot) * ...
                contact_coord(sample, y);   % y_C
            C(sample, 3) = C(sample, 3) + ( epsilon(k) / epsilon_tot) * ...
                contact_coord(sample, z);   % z_C
        end

        C(sample, 1) = alpha * C(sample, 1);
        C(sample, 2) = alpha * C(sample, 2);
        C(sample, 3) = (1 - alpha) * C(sample, 3);


        % z_G = 0.814;  
        % Accelerations
        if (sample > 1)
            x_G(sample, 3) = ...
                (1 / (M * ( Data(sample, 4) - C(sample, 3)))) * ...
                ( M * (ddZg + G) * (x_G(sample, 1) - ...
                C(sample, 1)) - dL(2) - tau_C_y );       % G x acceleration
            x_G(sample + 1, 2) = ...
                x_G(sample, 2) + x_G(sample, 3) * period;% G x speed
            x_G(sample + 1, 1) = ...
                2 * x_G(sample, 1) - x_G(sample - 1, 1) + ...
                period * period * x_G(sample, 3);        % G x position

            y_G(sample, 3) = ...
                (1 / (M * ( Data(sample, 4) - C(sample, 3)))) * ...             
                ( M * (ddZg + G) * (y_G(sample, 1) - ...
                C(sample, 2)) - dL(1) - tau_C_x );       % G y acceleration
            y_G(sample + 1, 2) = ...
                y_G(sample, 2) + y_G(sample, 3) * period;% G y speed
            y_G(sample + 1, 1) = ...
                2 * y_G(sample, 1) - y_G(sample - 1, 1) + ...
                period * period * y_G(sample, 3);        % G y position

        else
            x_G(1, 3) = 0;
            x_G(1, 2) = 0;
            x_G(1, 1) = Data(1, 2);

            x_G(2, 3) = 0;
            x_G(2, 2) = 0;
            x_G(2, 1) = Data(2, 2);

            y_G(1, 3) = 0;
            y_G(1, 2) = 0;
            y_G(1, 1) = Data(1, 3);

            y_G(2, 3) = 0;
            y_G(2, 2) = 0;
            y_G(2, 1) = Data(2, 3);        
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 7 : find  P_ref
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        P_ref_x = M * x_G(sample, 2);
        P_ref_y = M * y_G(sample, 2);
        P_ref_z = M * Data(sample, 8);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 8 : find xi_B_ref
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % xi_B = A^-1 y

        temp0 = zeros(6,6);
        temp0(1:3, 1:3) = M * I3;
        temp0(1:3, 4:6) = -M * hat(r_bc);
        temp0(4:6, 4:6) = I_tilde;



        temp1 = eye(6,6);
        temp1(1:3, 4:6) = -r_B_F1;

        temp2 = eye(6,6);
        temp2(1:3, 4:6) = -r_B_F2;

        temp3 = eye(6,6);
        temp3(1:3, 4:6) = -r_B_H1;

        temp4 = eye(6,6);
        temp4(1:3, 4:6) = -r_B_H2;


        A = temp0 ...
            - ([M_leg_1 ; H_leg_1] * J_leg_1\ temp1 + ...
            [M_leg_1 ; H_leg_1] * J_leg_2\ temp2) ...
            - ([M_leg_1 ; H_leg_1] * J_arm_1\ temp3 + ...
            [M_leg_1 ; H_leg_1] * J_arm_2\ temp4);

        y = [P_ref_x ; P_ref_y ; P_ref_z ; L] ...
            - ( [M_leg_1 ; H_leg_1] * J_leg_1\ xi_F_1 + ...
            [M_leg_1 ; H_leg_1] * J_leg_2\ xi_F_2 ) ...
            - ( [M_leg_1 ; H_leg_1] * J_arm_1\ xi_H_1 + ...
            [M_leg_1 ; H_leg_1] * J_arm_2\ xi_H_2 );

        xi_B_new = A\y;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 9 : checking convergence
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        if ( ( xi_B_new(1:2) - xi_B(1:2) ) < err )
            %fprintf('converge\n');
            converge = 1;
        else
            xi_ref = xi_B_new(1:2);
            xi = xi_B(1:2);
            %fprintf('ne converge pas\n');
            converge = 0;
        end


    end
    
    res_xi_B(sample, :) = xi_B_new;
    given_xi_B(sample, :) = xi_B;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 10 : find angular speeds with new linear and angular speed of B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % leg 1
    tmp = eye(6,6);
    tmp(1:3,4:6) = -r_B_F1;                % vector waist leg 1
    route = FindRoute(RLEG_JOINT5);
    J_leg_1 = CalcJacobian(route);
    d_theta_leg_1 = J_leg_1\ xi_F_1 - J_leg_1\ tmp * xi_B;  % angular speed

    % leg 2
    tmp = eye(6,6);
    tmp(1:3,4:6) = -r_B_F2;
    route = FindRoute(LLEG_JOINT5);
    J_leg_2 = CalcJacobian(route);
    d_theta_leg_2 = J_leg_2\ xi_F_2 - J_leg_2\ tmp * xi_B;

    % arm 1
    tmp = eye(6,6);
    tmp(1:3,4:6) = -r_B_H1;
    route = FindRoute(RARM_JOINT5);
    J_arm_1 = CalcJacobian(route(:,3:end));
    d_theta_arm_1 = J_arm_1\ xi_H_1 - J_arm_1\ tmp * xi_B;

    % arm 2
    tmp = eye(6,6);
    tmp(1:3,4:6) = -r_B_H2;
    route = FindRoute(LARM_JOINT5);
    J_arm_2 = CalcJacobian(route(:,3:end));
    d_theta_arm_2 = J_arm_2\ xi_H_2 - J_arm_2\ tmp * xi_B;

    dq(sample+1,:) = [d_theta_leg_1 ; d_theta_leg_2 ; zeros(4,1) ; d_theta_arm_1; 0 ; d_theta_arm_2 ; 0] ;
    q_values(sample+1,:) = q_values(sample,:) + dq(sample,:) * period;
    
    for i = 2:length(uLINK)
        uLINK(i).q   = q_values(i-1) ;
        uLINK(i).dq  = dq(i-1) ;
    end

    uLINK(WAIST).v = xi_B_new(1:3);
    uLINK(WAIST).w = xi_B_new(4:6);
    uLINK(WAIST).p = uLINK(WAIST).p + uLINK(WAIST).v * period ;
    
    
    
    ForwardKinematics(1);
    ForwardVelocity(1);
    
end



%hold on
%plot(res_xi_B(:, 1),res_xi_B(:, 2), 'r')
%plot(res_xi_B(:, 2), 'b')
%plot(x_G(:, 1), 'r')
%plot(y_G(:, 2), 'b')

%figure
%hold on
%plot(Whole_data(1:number_of_samples, 2), 'r')
%plot(Whole_data(1:number_of_samples, 3), 'b')
%plot(given_xi_B(:, 1), 'b')








