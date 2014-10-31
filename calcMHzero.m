
function [m_tilde, c_tilde, I_tilde] = calcMHzero(j)
global uLINK
  
  if j == 0
      m_tilde  = 0          ;
      c_tilde  = zeros(3,1) ;
      I_tilde  = zeros(3,3) ;
  else
      [m_tilde_sister , c_tilde_sister, I_tilde_sister] = calcMHzero(uLINK(j).sister);
      [m_tilde_child  , c_tilde_child , I_tilde_child ] = calcMHzero(uLINK(j).child );
  
      m_tilde  = uLINK(j).m + m_tilde_sister + m_tilde_child     ;                             % (23)  
      c_tilde  = (m_tilde_sister * c_tilde_sister + ...
                  m_tilde_child  * c_tilde_child  + ...
                  uLINK(j).m * (uLINK(j).p + uLINK(j).R * uLINK(j).c) ) / m_tilde ;
  
      I_tilde = I_tilde_sister + m_tilde_sister * D_(c_tilde_sister - c_tilde) + ...           % (25)
                I_tilde_child  + m_tilde_child  * D_(c_tilde_child  - c_tilde) + ...
                uLINK(j).R * uLINK(j).I * uLINK(j).R' + uLINK(j).m * D_(uLINK(j).c - c_tilde);
  
      mj = cross(uLINK(j).a , (c_tilde - uLINK(j).p))*m_tilde ;                                % (18)
      hj = cross(c_tilde , mj) + I_tilde * uLINK(j).a ;                                        % (19) 
  
      uLINK(j).m_tilde = m_tilde ;
      uLINK(j).c_tilde = c_tilde ;
      uLINK(j).I_tilde = I_tilde ;
      uLINK(j).mj      = mj      ;
      uLINK(j).hj      = hj      ;
      
      
    c1 = uLINK(j).R * uLINK(j).c;
    c  = uLINK(j).p + c1;
    P = uLINK(j).m * (uLINK(j).v + cross(c1, uLINK(j).w));
    L = cross(c, P) + uLINK(j).R * uLINK(j).I * uLINK(j).R' * uLINK(j).w;
    L = L + calcL(uLINK(j).sister) + calcL(uLINK(j).child);
      
  end
end

function ret = D_(r)
  r_temp = hat(r);
  ret = r_temp' * r_temp;
end