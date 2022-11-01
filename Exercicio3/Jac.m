function J = Jac(x)
  J(1,1) = 2 * x(1);
  J(1,2) = 2 * x(2);
  J(1,3) = 2 * x(3);
  J(2,1) = cos(x(1));
  J(2,2) = 3 * (x(2))^2;
  J(2,3) = 1;
  J(3,1) = 2 * x(1);
  J(3,2) = -2;
  J(3,3) = 1;
endfunction