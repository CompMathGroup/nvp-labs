#include "schemes.h"

#define R Rational

Alphas upstream(
	 R()     , R(1)    , R(),
R(), R(0, -1), R(-1, 1), R(), R(), 
R(), R()     , R()     , R(), R()
);

Alphas laxwendroff(
   	 R()           , R(1)       , R(),
R(), R(0, -.5, -.5), R(-1, 0, 1), R(0, .5, -.5), R(),
R(), R()           , R()        , R()          , R()
);

Alphas beamwarming(
   	           R()        , R(1)           , R(),
R(0, .5, -.5), R(0, -2, 1), R(-1, 1.5, -.5), R(), R(),
R()          , R()        , R()            , R(), R()
);

Alphas rusanov(
R(), R(1), R(),
R(0,1,0,-1,6,0,0,0),R(0,-1,-.5,.5),R(-1,.5,1,-.5),R(0,2,-3,1,6,0,0,0),R(),
R(), R(), R(), R(), R()
);

Alphas test(
R(0, -12, 6, 0, 12, -18, 0, 0),R(12, -6, -6, 0, 12, -18, 0, 0),R(),
R(0, 2, 3, 1, 12, -18, 0, 0),R(),R(-12, 12, 3, -3, 12, -18, 0, 0),
R(0, 4, -6, 2, 12, -18, 0, 0),R(),
R(),R(),R(),R(),R()
);

void selfTest() {
	bool all = true;
	std::cout << "Testing upstream" << std::endl;
	all &= upstream.sanity();
	std::cout << "Testing laxwendroff" << std::endl;
	all &= laxwendroff.sanity();
	std::cout << "Testing beamwarming" << std::endl;
	all &= beamwarming.sanity();
	std::cout << "Testing rusanov" << std::endl;
	all &= rusanov.sanity();
	std::cout << "Testing test" << std::endl;
	all &= test.sanity();
	std::cout << "All tests " << (all?"OK":"Failed") << std::endl;
}
