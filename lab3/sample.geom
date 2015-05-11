# Declare Robin conditions du/dn + 1 * u = 2
cond1 = robin { return 1; } { return 2; }
# And Neumann conditions du/dn = x * x + y
cond2 = neumann { return x * x + y; }

bigo = circle (0, 0) radius 5 neumann { return 0; }
eye1 = circle (-2, 1) radius 1 dirichlet { return 1; }
eye2 = circle (2, 1) radius 1 dirichlet { return -1; }

mouth = polygon (3, -1.5) dirichlet { reuturn x+y; } (-3, -1.5) cond2 (0, -2.5) cond2 close

nose1 = circle (1.5, 0.5) radius 1.8 dirichlet { return y; }
nose2 = circle (-1.5, 0.5) radius 1.8 dirichlet { return y; }

nose = nose1 intersect nose2
dom = bigo subtract (eye1 union eye2 union mouth union nose)

eq = elliptic { return uxx + uyy - 1; };

solve eq in dom meshsize 0.1 0.2 1
