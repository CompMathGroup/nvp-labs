%{
%}

%define api.namespace { config }
%define parser_class_name { Parser }
%locations
%initial-action {
	@$.begin.filename = @$.end.filename = &lexer.fname;
}
%code requires {
	#include <vector>
	#include "Value.h"
	#include "Reg.h"

	namespace config {

	typedef Register<Region> RegRegister;
	typedef Register<BoundaryCondition> BcRegister;
	typedef Register<Point> PointRegister;
	typedef Register<Problem> ProblemRegister;

	class Lexer;

	}
}
%code {
	#include "Lexer.h"

	using namespace std;

	static config::Parser::token_type
	yylex(config::Parser::semantic_type *yylval,
		config::Parser::location_type *yylloc,
		config::Lexer &lexer)
	{
		return lexer.lex(yylval, yylloc);
	}
}

%union {
	double dval;
	const Point *point;
	const std::string *sval;
	Path *path;
	const BoundaryCondition *bc;
	const Region *reg;
	const MeshSize *ms;
	const Config *ret;
	const Problem *prob;
}

%type<point> point cpoint
%type<path>  path
%type<bc>    bc cbc
%type<reg>	 region polygon circle cregion
%type<ms>	 meshsize
%type<prob>  problem cproblem
%type<ret>	 config

%token<dval> NUMBER
%token<sval> IDENTIFIER
%token<sval> CODE
%token       END 0
%token       POLYGON
%token       CLOSE
%token       CIRCLE
%token       RADIUS
%token       BNDCOND
%token       SOLVE
%token       ELLIPTIC
%token       MESHSIZE
%token       IN
%left        UNION SUBTRACT
%left        INTERSECT

%parse-param { config::Lexer &lexer}
%parse-param { const Config *&config }
%lex-param { config::Lexer &lexer }

%start config
%%

cproblem : ELLIPTIC CODE { $$ = new Problem(*$2); delete $2; ProblemRegister::add($$); }
		 ;

problem : cproblem   { $$ = $1; }
		| IDENTIFIER { $$ = ProblemRegister::lookup(*$1, @1); delete $1; }
		;

cpoint	: '(' NUMBER ',' NUMBER ')' { $$ = new Point($2, $4); PointRegister::add($$); }
		;

point	: cpoint
		| IDENTIFIER				{ $$ = PointRegister::lookup(*$1, @1); delete $1; }
		;

cbc		: BNDCOND CODE	{ $$ = new BoundaryCondition(*$2); delete $2; BcRegister::add($$); }
		;

bc		: cbc				{ $$ = $1; }
		| IDENTIFIER		{ $$ = BcRegister::lookup(*$1, @1); delete $1; }
		;

path	: path point bc		{ $1->append(*$2, $3); $$ = $1; }
		| point bc			{ $$ = new Path(); $$->append(*$1, $2); }
		;

polygon : POLYGON path CLOSE { $$ = new Polygon(*$2); delete $2; RegRegister::add($$); }
		;

circle	: CIRCLE point RADIUS NUMBER bc	{ $$ = new Circle(*$2, $4, $5); RegRegister::add($$); }
		;

cregion	: polygon					{ $$ = $1; }
		| circle					{ $$ = $1; }
		| region UNION region		{ $$ = new Union($1, $3); RegRegister::add($$); }
		| region SUBTRACT region	{ $$ = new Subtract($1, $3); RegRegister::add($$); }
		| region INTERSECT region	{ $$ = new Intersect($1, $3); RegRegister::add($$); }
		| '(' region ')'			{ $$ = $2; }
		;

region	: cregion					{ $$ = $1; }
		| IDENTIFIER				{ $$ = RegRegister::lookup(*$1, @1); delete $1; }
		;

statement : IDENTIFIER '=' cregion	{ RegRegister::add($3, *$1); delete $1; }
		  | IDENTIFIER '=' cbc		{ BcRegister::add($3, *$1); delete $1; }
		  | IDENTIFIER '=' cpoint	{ PointRegister::add($3, *$1); delete $1; }
		  | IDENTIFIER '=' cproblem	{ ProblemRegister::add($3, *$1); delete $1; }
		  ;

meshsize : MESHSIZE NUMBER NUMBER NUMBER { $$ = new MeshSize($2, $3, $4); }
		 | MESHSIZE NUMBER NUMBER        { $$ = new MeshSize($2, $3); }
		 | MESHSIZE NUMBER               { $$ = new MeshSize($2); }
		 ;

config	: SOLVE problem IN region meshsize { config = new Config($2, $4, *$5); delete $5; }
		| statement config	{ $$ = $2; }
		;
%%

void config::Parser::error(const config::location &lloc, const string &msg) {
	cerr << lloc << ": " << msg << endl;
}
