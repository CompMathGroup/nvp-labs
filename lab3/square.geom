dirichlet = bndcond { u }
neumann = bndcond { dudn }

dom = polygon (1, 1) dirichlet (-1, 1) dirichlet (-1, -1) dirichlet (1, -1) neumann close
eqn = elliptic { uxx + uyy - 1 }

solve eqn in dom meshsize 0.05
