prologues := 3;
%ahangle := 15;
%ahlength := 5;
outputtemplate := "%j-%c.eps";

beginfig(0);

a = 35;

b = 5;
z1 = (0,0);
z0 = z1 - (.5a,0);
z3 = z1 - (0,.5a);
z2 = z0 + (6a, 0);
z4 = z3 + (0, 3.5a);
z5 = z1 - (.5a,.5a);

picture axis;
def xtick(expr x, y) = draw (x - (0,.5b)) .. (x + (0,.5b)); label.bot(y, x - (0,.5b)); enddef;
def ytick(expr x, y) = draw (x - (.5b,0)) .. (x + (.5b,0)); label.lft(y, x - (.5b, 0)); enddef;
def dot(expr x) = filldraw fullcircle scaled 5 shifted x; enddef;
def odot(expr x) = draw fullcircle scaled 5 shifted x; enddef;
def diam(expr x) = filldraw ((5,0) -- (0,5) -- (-5,0) -- (0,-5) -- cycle) shifted x; enddef;

axis = image (
drawarrow z0 -- z2;
drawarrow z3 -- z4;

xtick (z0 + (a,0), btex $m-2$ etex);
xtick (z0 + (2a,0), btex $m-1$ etex);
xtick (z0 + (3a,0), btex $m$ etex);
xtick (z0 + (4a,0), btex $m+1$ etex);
xtick (z0 + (5a,0), btex $m+2$ etex);

ytick (z3 + (0,a), btex $n-1$ etex);
ytick (z3 + (0,2a), btex $n$ etex);
ytick (z3 + (0,3a), btex $n+1$ etex);
);

draw axis;

dot (z5 + (1a, a));
dot (z5 + (2a, a));
dot (z5 + (3a, a));
dot (z5 + (4a, a));
dot (z5 + (5a, a));

dot (z5 + (1a, 2a));
dot (z5 + (2a, 2a));
dot (z5 + (3a, 2a));
dot (z5 + (4a, 2a));
dot (z5 + (5a, 2a));

dot (z5 + (2a, 3a));
dot (z5 + (4a, 3a));

diam(z5 + (3a, 3a));

endfig;

beginfig(1);

draw axis;

z5 = -(.5a,.5a);

odot(z5 + (1a, a));
odot(z5 + (2a, a));
odot(z5 + (3a, a));
odot(z5 + (4a, a));
odot(z5 + (5a, a));

dot (z5 + (1a, 2a));
dot (z5 + (2a, 2a));
dot (z5 + (3a, 2a));
dot (z5 + (4a, 2a));
dot (z5 + (5a, 2a));

dot (z5 + (2a, 3a));
dot (z5 + (4a, 3a));

dot (z5 + (3a, 3a));

draw ((0,0) -- (a,0) -- (2a,0) -- (3a,0) -- (4a,0)) shifted (z5 + (a,2a));
draw ((2a,2a) -- (2a,3a)) shifted z5;
draw ((3a,2a) -- (3a,3a)) shifted z5;
draw ((4a,2a) -- (4a,3a)) shifted z5;

draw ((-5,-5) -- (5,5)) scaled .5 shifted (z5 + (3a,2.5a));
draw ((-5,5) -- (5,-5)) scaled .5 shifted (z5 + (3a,2.5a));

endfig;

beginfig(2);

draw axis;

z5 = -(.5a,.5a);

odot(z5 + (1a, a));
odot(z5 + (2a, a));
odot(z5 + (3a, a));
odot(z5 + (4a, a));
odot(z5 + (5a, a));

odot(z5 + (1a, 2a));
dot (z5 + (2a, 2a));
dot (z5 + (3a, 2a));
dot (z5 + (4a, 2a));
dot (z5 + (5a, 2a));

odot(z5 + (2a, 3a));
odot(z5 + (3a, 3a));
odot(z5 + (4a, 3a));


draw ((a,0) -- (2a,0) -- (3a,0) -- (4a,0)) shifted (z5 + (a,2a));

draw ((-5,-5) -- (5,5)) scaled .5 shifted (z5 + (3.5a,2a));
draw ((-5,5) -- (5,-5)) scaled .5 shifted (z5 + (3.5a,2a));

endfig;

end
