function err_norm = InverseKinmatics(to, Target)
global uLINK

lambda = 0.9;
ForwardKinematics(1);
idx = FindRoute(to);
for n = 1:10
  J   = CalcJacobian(idx);
  err = CalcVWerr(Target, uLINK(to));
  if norm(err) < 1E-6 break, end;
  dq = lambda * (J \ err);
  MoveJoints(idx, dq);
  ForwardKinematics(1);
end

if nargout == 1 
    err_norm = norm(err);
end