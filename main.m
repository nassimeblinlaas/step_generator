%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hirukawa multicontact step generator 2007
%
% Implémenté par Nassime BLIN 2014 pour Gepetto LAAS - CNRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc
clear







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delarations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global uLINK G
G = 9.81 ;



MG = 9.81 ;
pZ = 0.6487;            % position Z of robot constant

period = 0.005;         % sampling period in seconds

WAIST = 1;              % labels
RLEG_JOINT5 = 7;        %
LLEG_JOINT5 = 13;       %
RARM_JOINT5 = 23;       %
LARM_JOINT5 = 31;       %

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

fprintf('Reading ./morisawa.csv\n')
Data = csvread('./morisawa.csv');  
halfsitting = load('./halfsitting.dat');
fprintf('Reading done\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build contact points data
is_contact = zeros( length(Data) , 2);
contact_coord = zeros(length(Data), 4);
normal_vectors = zeros(length(Data), 3 * 2);    % three scalars for vector foreach two contact points

for i = 1:length(Data)
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


sample = sample + 1

M = TotalMass(1);	% total robot mass

uLINK(WAIST).p = ...
    [Data(sample + 1, 1), Data(sample + 1, 2), Data(sample + 1, 3)]' ;        % waist position
theta = Data(sample, 5);                                        % yaw = lacet = 0
uLINK(WAIST).R = [cos(theta),-sin(theta),0;...                  % waist rotation sin(theta),cos(theta),0;...
                  sin(theta),cos(theta),0;...
                  0,0,1 ];             
              
%uLINK(WAIST).v = ...
%    [Data(sample,6), Data(sample,7), Data(sample,8)]';         % waist speed
uLINK(WAIST).w = [0;0;0];                                       % waist angular speed

uLINK(RARM_JOINT5).p = uLINK(WAIST).p + [0.30;0;0];             % right arm position −30 cm
uLINK(LARM_JOINT5).p = uLINK(WAIST).p - [0.30;0;0];             % left  arm position +30 cm
uLINK(RLEG_JOINT5).p = ...
    [Data(sample,11), Data(sample,12), Data(sample,13)]';       % right foot position
uLINK(LLEG_JOINT5).p = ...
    [Data(sample,23), Data(sample,24), Data(sample,25)]';       % left foot position

for i = 1 : length(uLINK) - 1
    uLINK(i).q   = halfsitting(i);
end

ForwardKinematics(1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1 : give waist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_ref_B   = [ 0; 0; 0];     % waist speed vector
w_ref_B   = [ 0; 0; 0];
v_ref_F_1 = [ 0; 0; 0];
w_ref_H_1 = [ 0; 0; 0];
v_ref_F_2 = [ 0; 0; 0];
w_ref_H_2 = [ 0; 0; 0];


iteration   = 0;            % algorithm iteration
converge    = 0

while ( converge == 0 )

    iteration = iteration + 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2 : find angular speeds having new linear and angular speed of B
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % leg 1
    temp = eye(6,6);
    temp(1:3,4:6) = -r_B_F1;                % vector waist leg 1
    route = FindRoute(RLEG_JOINT5);
    J_leg_1 = CalcJacobian(route);
    J_leg_1_inv = inv(J_leg_1) * 200 / pi;
    d_theta_leg_1 = J_leg_1_inv * xi_F_1 - J_leg_1 * temp * xi_B;  % angular speeds

    % leg 2
    temp = eye(6,6);
    temp(1:3,4:6) = -r_B_F2;
    route = FindRoute(LLEG_JOINT5);
    J_leg_2 = CalcJacobian(route);
    J_leg_2_inv = inv(J_leg_2);
    d_theta_leg_2 = J_leg_2_inv * xi_F_2 - J_leg_2 * temp * xi_B;

    % bras 1
    temp = eye(6,6);
    temp(1:3,4:6) = -r_B_F1;
    route = FindRoute(RARM_JOINT5);
    J_arm_1 = CalcJacobian(route(:,3:end));
    J_arm_1_inv = inv(J_arm_1);
    d_theta_arm_1 = J_arm_1_inv * xi_H_1 - J_arm_1 * temp * xi_B;

    % bras 2
    temp = eye(6,6);
    temp(1:3,4:6) = -r_B_F1;
    route = FindRoute(LARM_JOINT5);
    J_arm_2 = CalcJacobian(route(:,4:end));
    J_arm_2_inv = inv(J_arm_2);
    d_theta_arm_2 = J_arm_2_inv * xi_H_2 - J_arm_2 * temp * xi_B;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [M_leg_1, H_leg_1, I_tilde_leg_1] = MH(RLEG_JOINT5);     %call MH subroutine
    [M_leg_2, H_leg_2, I_tilde_leg_2] = MH(LLEG_JOINT5);
    [M_arm_1, H_arm_1, I_tilde_arm_1] = MH(RARM_JOINT5);   
    [M_arm_2, H_arm_2, I_tilde_arm_2] = MH(LARM_JOINT5);

    d_theta = [d_theta_leg_1' d_theta_leg_2' d_theta_arm_1' d_theta_arm_2'];
    M_d_theta = [M_leg_1, M_leg_2, M_arm_1, M_arm_2];
    H_d_theta = [H_leg_1, H_leg_2, H_arm_1, H_arm_2];

    r_B_G = 0;              % vector waist to CoM


    I_tilde = I_tilde_arm_1 + I_tilde_arm_2 + I_tilde_leg_1 + I_tilde_leg_2;



    temp0 = zeros(6, 30);       
    temp0(1:3,1:3) = M*I3;
    temp0(1:3,4:6) = -M*r_B_G;
    temp0(4:6,4:6) = I_tilde;
    temp0(1:3,7:end) = M_d_theta;
    temp0(4:6,7:end) = H_d_theta;

    B = zeros(30,1);       
    B(1:3) = v_ref_B;
    B(4:6) = w_ref_B;
    B(7:end) = d_theta;

    PL = temp0 * B;
    P = PL(1:3);            % linear momentum
    L = PL(4:6);            % angular momentum



    dL = (L - L_prev) / period;     % finite difference method

    L_prev = L;             % save previous value for next step



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 4 : give \ddot{z_G}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % suppose that :
    ddZg = 0; % always



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5 : give lambda_k, slope of the ground
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % hypothesis : slope is zero
    % => Sum(n_kz)/K = 1
    % => alpha = 0   -> equations (16) (17) (18)
    alpha = 0;

    % or Somme(lambda_k) = 1

    % for two contact points :
    lambda = zeros(2, 1);
    lambda(1) = 0.5;
    lambda(2) = 0.5;

    epsilon = zeros(2, 1);
    epsilon_ref = zeros(2, 1);
    
    
    epsilon(1) = (1 - alpha) * M * ( ddZg + G ) * ...
        ( ( lambda(1) ) / ( lambda(1) * 1 + lambda(2) * 1 ) ) % multiplication by 1 is normal vectors n_k_z

    epsilon(2) = (1 - alpha) * M * ( ddZg + G ) * ...
        ( ( lambda(2) ) / ( lambda(2) * 1 + lambda(2) * 1 ) ) % multiplication by 1 is normal vectors n_k_z


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 6 : find d_d_x_G and d_d_y_G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    % Init
    tau_C_x = 0;
    tau_C_y = 0;
    epsilon_tot = 0;
    K = length(epsilon);
    C =   zeros(length(sample), 3);     % contact point x, y, z 
    x_G = zeros(length(sample), 3);     % for x_G, \dot x_G, \ddot x_G
    y_G = zeros(length(sample), 3);     % for y_G, \dot y_G, \ddot y_G
    
    
    % Torques
    for k = 1:K
        offset = (k-1) * 3;
        x = offset + 1;
        y = offset + 2;
        z = offset + 3;

        if ( is_contact(sample, k) == 1 )
            tau_C_x = tau_C_x + epsilon(k) * ( contact_coord(sample, y) * normal_vectors(sample, z) - ...
                contact_coord(sample, z) * normal_vectors(sample, y) );

            tau_C_y = tau_C_y + epsilon(k) * ( contact_coord(sample, x) * normal_vectors(sample, z) - ...
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
        
        C(sample, 1) = C(sample, 1) + ( epsilon(k) / epsilon_tot) * contact_coord(sample, x);   % x_C
        C(sample, 2) = C(sample, 2) + ( epsilon(k) / epsilon_tot) * contact_coord(sample, y);   % y_C
        C(sample, 3) = C(sample, 3) + ( epsilon(k) / epsilon_tot) * contact_coord(sample, z);   % z_C
    end
    
    C(sample, 1) = alpha * C(sample, 1);
    C(sample, 2) = alpha * C(sample, 2);
    C(sample, 3) = (1 - alpha) * C(sample, 3);
    
    
    % z_G = 0.814;  
    % Accelerations
    if (sample > 2)
        x_G(sample, 3) = (1 / (M * ( Data(sample, 4) - C(sample, 3)))) * ...
            ( M * (ddZg + G) * (x_G(sample, 1) - C(sample, 1)) - dL(2) - tau_C_y );      % G x acceleration
        x_G(sample + 1, 2) = x_G(sample, 2) + x_G(sample, 3) * period;                   % G x speed
        x_G(sample + 1, 1) = 2 * x_G(sample, 1) - x_G(sample - 1, 1) + period * period * x_G(sample, 3);    % G x position
            
        y_G(sample, 3) = (1 / (M * ( Data(sample, 4) - C(sample, 3)))) * ...             
            ( M * (ddZg + G) * (y_G(sample, 1) - C(sample, 2)) - dL(1) - tau_C_x );         % G x acceleration
        y_G(sample + 1, 2) = y_G(sample, 2) + y_G(sample, 3) * period;                      % G x speed
        y_G(sample + 1, 1) = 2 * y_G(sample, 1) - y_G(sample - 1, 1) + period * period * y_G(sample, 3);    % G x position
    else
        x_G(sample, 3) = 0
        x_G(sample, 2) = 0
        x_G(sample, 1) = Data(sample, 2);

        y_G(sample, 3) = 0
        y_G(sample, 2) = 0
        y_G(sample, 1) = Data(sample, 3);
    end
    
     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 7 : find  P_ref
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P_ref_x = M * x_G(sample, 2);
    P_ref_y = M * y_G(sample, 2);
    P_ref_z = M * Data(sample, 8);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 8 : find xi_B_ref
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % xi_B = A^-1 y

    M_F_1_star = [M_leg_1 ; H_leg_1] * J_leg_1_inv;
    M_F_2_star = [M_leg_1 ; H_leg_1] * J_leg_2_inv;
    M_H_1_star = [M_leg_1 ; H_leg_1] * J_arm_1_inv;
    M_H_2_star = [M_leg_1 ; H_leg_1] * J_arm_2_inv;


    temp0 = zeros(6,6);
    temp0(1:3, 1:3) = M * I3;
    temp0(1:3, 4:6) = -M * r_B_G;
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
        - (M_F_1_star * temp1 + M_F_2_star * temp2) ...
        - (M_H_1_star * temp3 + M_H_2_star * temp4);


    y = [P_ref_x ; P_ref_y ; P_ref_z ; L] ...
        - ( M_F_1_star * xi_F_1 + M_F_2_star * xi_F_2 ) ...
        - ( M_H_1_star * xi_H_1 + M_H_2_star * xi_H_2 );

    xi_B_ref = inv(A) * y;
    %xi_B_ref = A\y;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 9 : checking convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    if ( ( xi_B_ref(1:2) - xi_B(1:2) ) < err )
        fprintf('converge\n');
        converge = 1
    else
        fprintf('ne converge pas\n');
        converge = 0
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 10 : find angular speeds having new linear and angular speed of B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leg 1
temp = eye(6,6);
temp(1:3,4:6) = -r_B_F1;                % vector waist leg 1
route = FindRoute(RLEG_JOINT5);
J_leg_1 = CalcJacobian(route);
J_leg_1_inv = inv(J_leg_1) * 200 / pi;
d_theta_leg_1 = J_leg_1_inv * xi_F_1 - J_leg_1 * temp * xi_B_ref;  % angular speeds

% leg 2
temp = eye(6,6);
temp(1:3,4:6) = -r_B_F2;
route = FindRoute(LLEG_JOINT5);
J_leg_2 = CalcJacobian(route);
J_leg_2_inv = inv(J_leg_2);
d_theta_leg_2 = J_leg_2_inv * xi_F_2 - J_leg_2 * temp * xi_B_ref;

% bras 1
temp = eye(6,6);
temp(1:3,4:6) = -r_B_F1;
route = FindRoute(RARM_JOINT5);
J_arm_1 = CalcJacobian(route(:,3:end));
J_arm_1_inv = inv(J_arm_1);
d_theta_arm_1 = J_arm_1_inv * xi_H_1 - J_arm_1 * temp * xi_B_ref;

% bras 2
temp = eye(6,6);
temp(1:3,4:6) = -r_B_F1;
route = FindRoute(LARM_JOINT5);
J_arm_2 = CalcJacobian(route(:,4:end));
J_arm_2_inv = inv(J_arm_2);
d_theta_arm_2 = J_arm_2_inv * xi_H_2 - J_arm_2 * temp * xi_B_ref;






