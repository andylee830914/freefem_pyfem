#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <sstream>
#include <dlfcn.h>
#include <vector>
#ifndef WITH_NO_INIT
#include "ff++.hpp"
#include "AFunction.hpp"
#include "AFunction_ext.hpp"
#endif
using namespace std;
using namespace Fem2D;
namespace py = pybind11;
using namespace py::literals;

namespace{
    const char DEFAULT_FETYPE[] = "P1";
    const char DEFAULT_CMM[] = "";
    const char DEFAULT_SCRIPT[] = "";
    const char BASE_DIR[] = ".";
    int plotcount = 0;
    void *python;
    vector<vector<int>> tri;
    vector<vector<double>> mesh;
}



std::string get_string(Stack stack, Expression e, const char *const DEFAULT)
{
    const size_t length = 128;
    char *const carg = new char[length];

    if (!e)
        strcpy(carg, DEFAULT);
    else
        strncpy(carg, GetAny<string *>((*e)(stack))->c_str(), length);

    return std::string(carg);
}


class WEBPLOT_Op : public E_F0mps
{
  public:
    Expression eTh, ef, empi;
    static const int n_name_param = 3;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];

    double arg(int i, Stack stack, double defvalue) const { return nargs[i] ? GetAny<double>((*nargs[i])(stack)) : defvalue; }
    long arg(int i, Stack stack, long defvalue) const { return nargs[i] ? GetAny<long>((*nargs[i])(stack)) : defvalue; }
    KN<double> *arg(int i, Stack stack, KN<double> *defvalue) const { return nargs[i] ? GetAny<KN<double> *>((*nargs[i])(stack)) : defvalue; }
    bool arg(int i, Stack stack, bool defvalue) const { return nargs[i] ? GetAny<bool>((*nargs[i])(stack)) : defvalue; }

  public:
    WEBPLOT_Op(const basicAC_F0 &args, Expression th)
    : eTh(th), ef(0)
    {
        args.SetNameParam(n_name_param, name_param, nargs);
    }

    WEBPLOT_Op(const basicAC_F0 &args, Expression f, Expression th)
    : eTh(th), ef(f)
    {
        args.SetNameParam(n_name_param, name_param, nargs);
    }

    AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type WEBPLOT_Op::name_param[] =
    {
        // modify static const int n_name_param = ... in the above member
        {"cmm", &typeid(string *)},
        {"fetype", &typeid(string *)},
        {"script", &typeid(string *)}
        //{  "logscale",  &typeid(bool)} // not implemented
};

AnyType WEBPLOT_Op::operator()(Stack stack) const
{
    const std::string cmm = get_string(stack, nargs[0], DEFAULT_CMM);
    const std::string fetype = get_string(stack, nargs[1], DEFAULT_FETYPE);
    const std::string script = get_string(stack, nargs[2], DEFAULT_SCRIPT);

    const Mesh *const pTh = GetAny<const Mesh *const>((*eTh)(stack));
    ffassert(pTh);
    const Fem2D::Mesh &Th(*pTh);
    const int nVertices = Th.nv;
    const int nTriangles = Th.nt;

    R2 Pmin, Pmax;
    Th.BoundingBox(Pmin, Pmax);

    const double &x0 = Pmin.x;
    const double &y0 = Pmin.y;

    const double &x1 = Pmax.x;
    const double &y1 = Pmax.y;


    const double unset = -1e300;


    KN<double> f_FE(Th.nv, unset);

    if (fetype != "P1")
    {
        std::cout << "plotPDF() : Unknown fetype : " << fetype << std::endl;
        std::cout << "plotPDF() : Interpolated as P1 (piecewise-linear)" << std::endl;
    }
    mesh.resize(3 * Th.nt);
    for (int it = 0; it < Th.nt; it++)
    {
        vector<int> v;
        for (int iv = 0; iv < 3; iv++)
        {
            vector<double> m;
            int i = Th(it, iv);
            int j = iv + it * 3;
            // std::cout << j << std::endl;
            
            v.push_back(j);

            MeshPointStack(stack)->setP(pTh, it, iv);
            double temp;
            if (ef)
            {
                temp = GetAny<double>((*ef)(stack)); // Expression ef is atype<double>()
            }
            else
            {
                temp = 0;
            }

            m.push_back(Th(i).x);
            m.push_back(Th(i).y);
            m.push_back(temp);

            mesh[j]= m;
        }
        tri.push_back(v);
    }


    python = dlopen("/home/andylee/.pyenv/versions/3.10.4/lib/libpython3.10.so", RTLD_NOW | RTLD_GLOBAL);

    py::scoped_interpreter guard{};
    auto py_module = py::module_::import(script.c_str());

    py::finalize_interpreter();
    dlclose(python);
    return true;
}


class WEBPLOT : public OneOperator
{
    const int argc;

  public:
    WEBPLOT()    : OneOperator(atype<long>(),                  atype<const Mesh *>()), argc(1) {}
    WEBPLOT(int) : OneOperator(atype<long>(), atype<double>(), atype<const Mesh *>()), argc(2) {}

    E_F0 *code(const basicAC_F0 &args) const
    {
        if (argc == 1)
            return new WEBPLOT_Op(args, t[0]->CastTo(args[0]));
        else if (argc == 2)
            return new WEBPLOT_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        else
            ffassert(0);
    }
};

// convert a vector of vectors to a numpy array
template <typename T>
py::array_t<T> convert_to_numpy(const std::vector<std::vector<T>> &data)
{
    size_t rows = data.size();
    size_t cols = (rows > 0) ? data[0].size() : 0;

    py::array_t<T> result = py::array_t<T>({rows, cols});

    auto buffer_info = result.request();
    T *ptr = static_cast<T *>(buffer_info.ptr);

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            ptr[i * cols + j] = data[i][j];
        }
    }

    return result;
}

PYBIND11_EMBEDDED_MODULE(freefem, m)
{
    m.attr("fetype") = DEFAULT_FETYPE;
    m.attr("tri") = convert_to_numpy(tri);
    m.attr("mesh") = convert_to_numpy(mesh);
}

static void init(){
    Global.Add("python", "(", new WEBPLOT());
    Global.Add("python", "(", new WEBPLOT(0));
}

LOADFUNC(init);
