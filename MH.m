%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcule les matrices d'inertie linéraire et angulaire pour chaque joint à
% partir de l'extrémité fournie en argument.
%
% Algorithme issu de "Resolved Momentum Control: Humanoid Motion Planning
% based on the Linear and Angular Momentum" par KAJITA 2003.
%
% Implémenté par Nassime BLIN 2014 pour Gepetto LAAS - CNRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [M_theta, H_theta, I_tilde] = MH(joint)
    global uLINK

    %clc
    
    m_j = 0;
    m_j_1 = 0;
    m_tilde = 0;
    c_tilde = [0;0;0];
    I_tilde = 0;
    m_tilde_j_1 = 0;
    c_tilde_j_1 = 0;
    I_tilde_j_1 = 0;
    
    
    route = FindRoute(joint);
    if length(route)>6
        %disp('changement')
        route = route(length(route)-5:end);
    end
    
    i = 0;
    for j = joint : -1 : route(1)
        
        i = i + 1; 
        
        m_j_1 = uLINK(j).m;
        c_j_1 = uLINK(j).c;
        R_j_1 = uLINK(j).R;
        I_j_1 = uLINK(j).I;
        a_j = uLINK(j).a;
        r_j = uLINK(j).p;
        
        m_tilde_j_1 = m_tilde + m_j_1;                  % (23)
        c_tilde_j_1 = (m_tilde * c_tilde + m_j_1 * c_j_1) / (m_tilde + m_j_1);  % (24)
        I_tilde_j_1 = ...
            I_tilde + ...
            m_tilde * D(c_j_1 - c_tilde_j_1) + ...
            R_j_1 * I_j_1 * R_j_1' + ...
            m_j_1 * D(c_j_1 - c_tilde_j_1);             % (25)




        m_j = cross(a_j , (c_tilde - r_j) * m_tilde);   % (18)
        h_j = cross(c_tilde , m_j) + I_tilde * a_j;     % (19)
        
        M_theta(:,i) = m_j;
        H_theta_zero(:,i) = h_j;
        
        
        % save previous result
        m_tilde = m_tilde_j_1;
        c_tilde = c_tilde_j_1;
        I_tilde = I_tilde_j_1;   
        
        
    end

    H_theta = H_theta_zero - hat(c_tilde) * M_theta;
    
end


function ret = D(r)

    r_temp = hat(r);
    ret = r_temp' * r_temp;
end