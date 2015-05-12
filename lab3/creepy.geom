cond1 = bndcond { dudn + u - 2 }
cond2 = bndcond { dudn - x * x - y }

bigo = circle (0, 0) radius 5 bndcond { dudn }
eye1 = circle (-2, 1) radius 1 bndcond { u + 1 }
eye2 = circle (2, 1) radius 1 bndcond { u - 1 }

mouth = polygon (3, -1.5) cond1 (-3, -1.5) cond2 (0, -2.5) cond2 close

nose1 = circle (1.5, 0.5) radius 1.8 bndcond { y - u }
nose2 = circle (-1.5, 0.5) radius 1.8 bndcond { y - u }

nose = nose1 intersect nose2
dom = bigo subtract (eye1 union eye2 union mouth union nose)

eq = elliptic { uxx + uyy }

solve eq in dom meshsize 0.4 0.2 0.6
