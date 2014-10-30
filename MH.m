%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcule les matrices d'inertie linéraire et angulaire pour chaque joint à
% partir de la chaîne cinématique fournie en argument.
%
% Algorithme issu de "Resolved Momentum Control: Humanoid Motion Planning
% based on the Linear and Angular Momentum" par KAJITA 2003.
%
% Implémenté par Nassime BLIN 2014 pour Gepetto LAAS - CNRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ M_d_theta, H_d_theta, I_tilde_waist, M_leg_1, M_leg_2, ...
            M_arm_1, M_arm_2, H_leg_1, H_leg_2, H_arm_1, H_arm_2 ] = MH()
    global uLINK
    
    WAIST = 1;         %labels
    RLEG_JOINT0 = 2;   %
    RLEG_JOINT1 = 3;   %
    RLEG_JOINT2 = 4;   %
    RLEG_JOINT3 = 5;   %
    RLEG_JOINT4 = 6;   %
    RLEG_JOINT5 = 7;   %
    LLEG_JOINT0 = 8;   %
    LLEG_JOINT1 = 9;   %
    LLEG_JOINT2 = 10;  %
    LLEG_JOINT3 = 11;  %
    LLEG_JOINT4 = 12;  %
    LLEG_JOINT5 = 13;  %
    CHEST_JOINT0 = 14; %
    CHEST_JOINT1 = 15; %
    HEAD_JOINT0 = 16;  %
    HEAD_JOINT1 = 17;  %
    RARM_JOINT0 = 18;  %
    RARM_JOINT1 = 19;  %
    RARM_JOINT2 = 20;  %
    RARM_JOINT3 = 21;  %
    RARM_JOINT4 = 22;  %
    RARM_JOINT5 = 23;  %
    RARM_JOINT6 = 24;  %
    LARM_JOINT0 = 25;  %
    LARM_JOINT1 = 26;  %
    LARM_JOINT2 = 27;  %
    LARM_JOINT3 = 28;  %
    LARM_JOINT4 = 29;  %
    LARM_JOINT5 = 30;  %
    LARM_JOINT6 = 31;  %

    route = FindRoute(HEAD_JOINT1);    % head
    route_head = route(3:4);
    route_leg_1 = FindRoute(RLEG_JOINT5);    % Rleg
    route_leg_2 = FindRoute(RLEG_JOINT5);    % Lleg
    route = FindRoute(RARM_JOINT5);    % Rarm
    route_arm_1 = route(3:end);
    route = FindRoute(LARM_JOINT5);    % Larm
    route_arm_2 = route(3:end);
    
    [M_leg_1, H_leg_1, I_tilde_leg_1, m_tilde_leg_1, c_tilde_leg_1] = MH_limb(route_leg_1);     % call MH subroutine
    [M_leg_2, H_leg_2, I_tilde_leg_2, m_tilde_leg_2, c_tilde_leg_2] = MH_limb(route_leg_2);
    [M_arm_1, H_arm_1, I_tilde_arm_1, m_tilde_arm_1, c_tilde_arm_1] = MH_limb(route_arm_1);   
    [M_arm_2, H_arm_2, I_tilde_arm_2, m_tilde_arm_2, c_tilde_arm_2] = MH_limb(route_arm_2);
    [M_head, H_head, I_tilde_head, m_tilde_head, c_tilde_head] = MH_limb(route_head);

    
   
    % Chest Inertia Matrices
    m_tilde_chest_1 = m_tilde_arm_1 + m_tilde_arm_2 + m_tilde_head + uLINK(CHEST_JOINT1).m;
    c_tilde_chest_1 = (m_tilde_arm_1*c_tilde_arm_1 + m_tilde_arm_2*c_tilde_arm_2 + m_tilde_head*c_tilde_head + ...
        uLINK(CHEST_JOINT1).m * uLINK(CHEST_JOINT1).c) / m_tilde_chest_1 ;
    m_tilde_chest_0 = m_tilde_chest_1 + uLINK(CHEST_JOINT0).m;
    c_tilde_chest_0 = (m_tilde_chest_1*c_tilde_chest_1 + ...
        uLINK(CHEST_JOINT0).m * uLINK(CHEST_JOINT0).c) / m_tilde_chest_0 ;
    
    I_tilde_chest_1 = ...
        I_tilde_arm_1 + m_tilde_arm_1 * D(c_tilde_arm_1 - c_tilde_chest_1) + ...
        I_tilde_arm_2 + m_tilde_arm_2 * D(c_tilde_arm_2 - c_tilde_chest_1) + ...
        I_tilde_head + m_tilde_head   * D(c_tilde_head -  c_tilde_chest_1) + ...
        uLINK(CHEST_JOINT1).R * uLINK(CHEST_JOINT1).I * uLINK(CHEST_JOINT1).R' + ...
            uLINK(CHEST_JOINT1).m * D(uLINK(CHEST_JOINT1).c - c_tilde_chest_1);

    I_tilde_chest_0 = ...
        I_tilde_chest_1 + ...
        m_tilde_chest_1 * D(c_tilde_chest_1 - c_tilde_chest_0) + ...
        uLINK(CHEST_JOINT0).R * uLINK(CHEST_JOINT0).I * uLINK(CHEST_JOINT0).R' + ...
        uLINK(CHEST_JOINT0).m * D(uLINK(CHEST_JOINT0).c - c_tilde_chest_0);
    
    M_chest(:, 2) = cross(uLINK(CHEST_JOINT1).a , (c_tilde_chest_1 - uLINK(CHEST_JOINT1).p))*m_tilde_chest_1;
    M_chest(:, 1) = cross(uLINK(CHEST_JOINT0).a , (c_tilde_chest_0 - uLINK(CHEST_JOINT0).p))*m_tilde_chest_0;
    H_zero_chest(:, 2) = cross(c_tilde_chest_1 , M_chest(:, 2)) + I_tilde_chest_1 * uLINK(CHEST_JOINT1).a;      % (19)
    H_zero_chest(:, 1) = cross(c_tilde_chest_0 , M_chest(:, 1)) + I_tilde_chest_0 * uLINK(CHEST_JOINT0).a;     % (19)
    H_chest = H_zero_chest - hat(I_tilde_chest_0) * M_chest;
    
    
    
    % Waist Inertia Matrices
    m_tilde_waist = m_tilde_leg_1 + m_tilde_leg_2 + m_tilde_chest_0 + uLINK(WAIST).m;
    c_tilde_waist = (m_tilde_leg_1*c_tilde_leg_1 + m_tilde_leg_2*c_tilde_leg_2 + m_tilde_chest_0*c_tilde_chest_0 + ...
        uLINK(WAIST).m * uLINK(WAIST).c) / m_tilde_waist ;

    I_tilde_waist = ...
        I_tilde_leg_1   + m_tilde_leg_1   * D(c_tilde_leg_1 - c_tilde_waist) + ...
        I_tilde_leg_2   + m_tilde_leg_2   * D(c_tilde_leg_2 - c_tilde_waist) + ...
        I_tilde_chest_0 + m_tilde_chest_0 * D(c_tilde_chest_0 - c_tilde_waist) + ...
        uLINK(WAIST).R * uLINK(WAIST).I * uLINK(WAIST).R' + ...
            uLINK(WAIST).m * D(uLINK(WAIST).c - c_tilde_waist);
    
    M_waist = cross(uLINK(WAIST).a , (c_tilde_waist - uLINK(WAIST).p))*m_tilde_waist ;
    H_zero_waist = cross(c_tilde_waist , M_waist) + I_tilde_waist * uLINK(WAIST).a ;     % (19)
    H_waist = H_zero_chest - hat(I_tilde_waist) * M_chest;
    
    
    
    % put everything together
    M_d_theta = [M_leg_1, M_leg_2, M_arm_1, M_arm_2];
    H_d_theta = [H_leg_1, H_leg_2, H_arm_1, H_arm_2];


end


function [M_theta, H_theta, I_tilde_j, m_tilde_j, c_tilde_j] = MH_limb(route)
    global uLINK
    
    
    
    m_tilde_j = uLINK(route(end)).m;
    c_tilde_j = uLINK(route(end)).c;
    I_tilde_j = uLINK(route(end)).I; 
    a_j     = uLINK(route(end)).a;
    r_j     = uLINK(route(end)).p;
    
    length(route);
    
    M_theta(:, 1) = cross(a_j , (c_tilde_j - r_j) * m_tilde_j);             % (18)
    H_theta_zero(:, 1) = cross(c_tilde_j , M_theta(:, 1)) + I_tilde_j * a_j;     % (19)
    
    for ii = length(route): - 1 : 2
        
        j = route(ii);
        j_1 = route(ii-1);
            
        m_j_1   = uLINK(j_1).m;
        c_j_1   = uLINK(j_1).c;
        R_j_1   = uLINK(j_1).R;
        I_j_1   = uLINK(j_1).I;
        a_j     = uLINK(j).a;
        r_j     = uLINK(j).p;
  
        
        m_tilde_j_1 = m_tilde_j + m_j_1;                                              % (23)
        c_tilde_j_1 = (m_tilde_j * c_tilde_j + m_j_1 * c_j_1) / (m_tilde_j + m_j_1);  % (24)
        I_tilde_j_1 = ...
            I_tilde_j + ...
            m_tilde_j * D(c_tilde_j - c_tilde_j_1) + ...
            R_j_1 * I_j_1 * R_j_1' + ...
            m_j_1 * D(c_j_1 - c_tilde_j_1);                                           % (25)


        m_j = cross(a_j , (c_tilde_j - r_j) * m_tilde_j);   % (18)
        h_j = cross(c_tilde_j , m_j) + I_tilde_j * a_j;     % (19)
        
        M_theta(:, ii) = m_j;
        H_theta_zero(:, ii) = h_j;
        
        
        % save previous result
        m_tilde_j = m_tilde_j_1;
        c_tilde_j = c_tilde_j_1;
        I_tilde_j = I_tilde_j_1;   
        
        
    end

    
    H_theta = H_theta_zero - hat(c_tilde_j) * M_theta;
    
end


function ret = D(r)

    r_temp = hat(r);
    ret = r_temp' * r_temp;
end

