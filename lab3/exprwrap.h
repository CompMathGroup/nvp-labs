#ifndef __EXPR_WRAP_H__
#define __EXPR_WRAP_H__

#include <string>
#include <memory>
#include <cassert>
#include <vector>

namespace interpreted {

struct ExpressionImpl;

class Expression {
    std::unique_ptr<ExpressionImpl> pimpl;
    mutable std::vector<double> args;
    double value() const;
    template<int m>
    double call(double arg) const
    {
        assert(m == args.size() - 1);
        args[m] = arg;
        return value();
    }
    template<int m, typename ...Args>
    double call(double arg, Args ...args) const {
        this->args[m] = arg;
        return call<m + 1>(args...);
    }
public:
    Expression(const std::string &s, const std::vector<std::string> &argnames);
    ~Expression();
    template<typename ...Args>
    double operator()(Args ...args) const {
        return call<0>(args...);
    }
};

}

#endif
