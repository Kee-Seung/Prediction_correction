function alpha = Backtracking_Line_Search(obj,f0val,x,df0dx,Dobj_old,dir,dir_old,alpha_old)


alpha_ini = alpha_old*Dobj_old'*dir_old/(df0dx'*dir); % Initial step length
alpha = alpha_ini;
fac = 1/2; % reduction factor of alpha: typically 0.5 is used.
c = 1e-4; % Typical value

while obj(x+alpha*dir) > f0val + c*alpha*dir'*df0dx

    alpha = fac*alpha;

end
