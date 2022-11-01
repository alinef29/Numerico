function [x,res] = nonlinNewton(fun ,Jac, x0, tol, nmax)
niter = 0;
err = tol + 1;
x = x0;

while err >= tol & niter < nmax
  J = Jac(x);
  F = fun(x);
  deltax = - J\F;
  x = x + deltax;
  err = norm(deltax);
  niter = niter + 1;
end
res = norm(fun(x));
if (niter== nmax & err > tol)
  fprintf ([' Fails to converge within maximum number of iterations .\n The iterate returned has relative residual %e\n'],F);
else
  fprintf (['The method converged at iteration %i with residual %e\n'],niter ,F);
end
return
