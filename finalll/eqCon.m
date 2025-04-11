function [c, ceq] = eqCon(x)

ceq = [];
rad = 8;
tol = 1e-3;
confcnval = sum(x) - rad;
c = [confcnval - tol;-confcnval - tol];