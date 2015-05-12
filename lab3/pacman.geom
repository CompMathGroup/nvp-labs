neumann = bndcond { dudn }
dirichlet = bndcond { u }

o = circle (0, 0) radius 1 neumann
cut = polygon (2, 1) dirichlet (0, 0) dirichlet (2, -1) dirichlet close
eye = circle (0, 0.6) radius 0.2 neumann
dom = o - cut - eye

eqn = elliptic { uxx + uyy - 1 }

solve eqn in dom meshsize 0.07 0.03 1
