outer = circle (0, 0) radius 4 bndcond { u - 4 * sin(5 * atan2(y, x)) }
inner = circle (0, 0) radius 2 bndcond { u }

solve elliptic { uxx + uyy } in (outer - inner) meshsize 0.2
