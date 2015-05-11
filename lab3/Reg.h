#ifndef __REG_H__
#define __REG_H__

#include <map>
#include <list>
#include <memory>
#include <string>
#include <iostream>
#include <sstream>

#include "location.hh"

using loc = config::location;

template<class T>
class Register {
    typedef std::map<std::string, const T *> Named;
    typedef std::list<const T *> Unnamed;

    static void purge_named(Named *v) {
        // Named objects once were unnamed, so are deleted in purge_unnamed
        delete v;
    }

    static void purge_unnamed(Unnamed *v) {
        for (const auto &it : *v) {
            delete it;
        }
        delete v;
    }

    static Named &named() {
        static std::unique_ptr<Named, void(*)(Named *)> r(new Named(), purge_named);
        return *r.get();
    }
    static Unnamed &unnamed() {
        static std::unique_ptr<Unnamed, void(*)(Unnamed *)> r(new Unnamed(), purge_unnamed);
        return *r.get();
    }
public:
    static const T *lookup(const std::string &name, const loc &lloc) {
        const auto it = named().find(name);
        if (it == named().end()) {
            std::stringstream ss;
            ss << lloc;
            std::cerr << ss.str() << ": Could not find `" << name << "'" << std::endl;
            throw std::out_of_range("No `" + name + "' in " + __PRETTY_FUNCTION__);
        }
        return it->second;
    }
    static void add(const T *p, const std::string name = "") {
//        std::cout << "Adding `" << name << "' p = " << p << std::endl;
        if (name == "")
            unnamed().push_back(p);
        else
            named()[name] = p;
    }
};

#endif
