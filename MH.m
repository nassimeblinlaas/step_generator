%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcule les matrices d'inertie linéraire et angulaire pour chaque joint à
% partir de l'extrémité fournie en argument.
%
% Algorithme issu de "Resolved Momentum Control: Humanoid Motion Planning
% based on the Linear and Angular Momentum" par KAJITA 2003.
%
% Implémenté par Nassime BLIN 2014 pour Gepetto LAAS - CNRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [M_theta, H_theta, I_tilde_j] = MH(joint)
    global uLINK

    route = FindRoute(joint);
    
    m_tilde_j = uLINK(route(length(route))).m;
    c_tilde_j = uLINK(route(length(route))).c;
    I_tilde_j = uLINK(route(length(route))).I;    
    
    for ii = length(route): - 1 : 1
        
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
        
        M_theta(:,i) = m_j;
        H_theta_zero(:,i) = h_j;
        
        
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