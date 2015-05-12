#include "exprwrap.h"

#include <exprtk.hpp>

#include <vector>

namespace interpreted {

using symbol_table = exprtk::symbol_table<double>;
using expression = exprtk::expression<double>;
using parser = exprtk::parser<double>;

struct ExpressionImpl {
    symbol_table vars;
    expression expr;
};

Expression::Expression(const std::string &s, const std::vector<std::string> &argnames)
    : pimpl(new ExpressionImpl), args(argnames.size())
{
    for (size_t i = 0; i < args.size(); i++) {
        pimpl->vars.add_variable(argnames[i], args[i]);
    }
    pimpl->vars.add_constants();
    pimpl->expr.register_symbol_table(pimpl->vars);

    parser p;
    if (!p.compile(s, pimpl->expr))
        throw std::invalid_argument("Failed to compile expression `" + s + "'");
}

Expression::~Expression() {
}

double Expression::value() const {
    return pimpl->expr.value();
}

}
