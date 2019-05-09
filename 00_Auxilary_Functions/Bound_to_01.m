% This function bouds the value in X between 0 and 1 while preserving their
% ranks. If Tol = .01 then it keeps everything between for example .01 and .99 as they are and uses a
% nonlinear transformation to squeez in the rest of the values between
% 0-.01 or .99-1

function Y = Bound_to_01(X)

Tol = .01; % this is the length of the area that will be used for nonlinear transformation

Lb = Tol; % lower bound
Ub = 1-Tol; % upper bound

Fl = X<Lb;
Fu = X>Ub;


Y = X;
Y(Fl) = Lb + Tol*transform(X(Fl)-Tol);
Y(Fu) = Ub + Tol*transform(X(Fu)-Tol);


end



function t = transform(d)
t = Pej_Transform_logistic(d)-.5;
end