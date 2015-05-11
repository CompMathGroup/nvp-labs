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
	Region *ret;
}

%type<point> point cpoint
%type<path>  path
%type<bc>    bc cbc
%type<reg>	 region polygon circle cregion config

%token<dval> NUMBER
%token<sval> IDENTIFIER
%token<sval> CODE
%token       END 0
%token       POLYGON
%token       CLOSE
%token       CIRCLE
%token       RADIUS
%token       NEUMANN
%token       DIRICHLET
%token       ROBIN
%left        UNION SUBTRACT
%left        INTERSECT

%parse-param { config::Lexer &lexer}
%parse-param { const Region *&config }
%lex-param { config::Lexer &lexer }

%start config
%%

cpoint	: '(' NUMBER ',' NUMBER ')' { $$ = new Point($2, $4); PointRegister::add($$); }
		;

point	: cpoint
		| IDENTIFIER				{ $$ = PointRegister::lookup(*$1, @1); delete $1; }
		;

cbc		: DIRICHLET CODE	{ $$ = new Diriclet(*$2); delete $2; BcRegister::add($$); }
		| NEUMANN CODE		{ $$ = new Neumann(*$2); delete $2; BcRegister::add($$); }
		| ROBIN CODE CODE	{ $$ = new Robin(*$2, *$3); delete $2; delete $3; BcRegister::add($$); }
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

statement : IDENTIFIER '=' cregion	{ RegRegister::add($3, *$1); delete $1; }
		  | IDENTIFIER '=' cbc		{ BcRegister::add($3, *$1); delete $1; }
		  | IDENTIFIER '=' cpoint	{ PointRegister::add($3, *$1); delete $1; }
		  ;

config	: region			{ config = $1; }
		| statement config	{ $$ = $2; }
		;
%%

void config::Parser::error(const config::location &lloc, const string &msg) {
	cerr << lloc << ": " << msg << endl;
}
