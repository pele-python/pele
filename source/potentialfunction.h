#ifndef _PELE_POTENTIAL_FUNCTION_H
#define _PELE_POTENTIAL_FUNCTION_H
#include "array.h"
#include "base_potential.h"
#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <stdexcept>

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

using std::cout;
using std::cerr;

namespace pele {

    /**
     * This class wraps a get_energy function and a get_energy_gradient function in
     * a class that derives from BasePotential.  This is necessary to be able to use
     * the functions in the pele c++ interface.  This is a backup method, the
     * preferred method is to define a class separately for each potential, which
     * eliminates the need to pass around the void * userdata parameter.
     */
    class PotentialFunction : public BasePotential
    {
        public:
            typedef double EnergyCallback(Array<double>, void *);
            typedef double EnergyGradientCallback(Array<double>, Array<double>, void *);

            PotentialFunction(EnergyCallback *get_energy, EnergyGradientCallback *get_energy_gradient, void *userdata)
                :    _get_energy(get_energy), _get_energy_gradient(get_energy_gradient), _userdata(userdata) {}

            virtual double get_energy(Array<double> x) { return (*_get_energy)(x, _userdata); } ;
            virtual double get_energy_gradient(Array<double> x, Array<double> grad) {  return (*_get_energy_gradient)(x, grad, _userdata); }

        private:
            EnergyCallback *_get_energy;
            EnergyGradientCallback *_get_energy_gradient;
            void *_userdata;
    };


    class PythonPotential : public BasePotential
    {
            PyObject * _potential;

        public:
            PythonPotential(PyObject * potential)
                :    _potential(potential)
            {
                Py_XINCREF(_potential);
                import_array();
                std::cout << "in constructor: done calling numpy import_array()\n";
            }
            ~PythonPotential() { 
                Py_XDECREF(_potential); 
            }

            virtual double get_energy(Array<double> x) 
            { 
                // create a numpy array from x
                // this wraps the data which is unsafe becase the python object might live longer than the data in x.data
                std::cout << "creating numpy array\n";
                npy_intp N = (npy_intp) x.size();
                PyObject * numpyx = PyArray_SimpleNewFromData(1, &N, NPY_DOUBLE, static_cast<void*>(x.data()));
                if (!numpyx){
                    cerr << "created numpy object is NULL\n";
                    throw std::runtime_error("created numpy object is NULL\n");
                }
                
                // call the function getEnergy
                PyObject * name = PyString_FromString("getEnergy");
                PyObject * returnval = PyObject_CallMethodObjArgs(_potential, name, numpyx, NULL);
                Py_XDECREF(name); 
                Py_XDECREF(numpyx); 
                if (!returnval){
                    //parse error
                    cout << "getEnergy returned a NULL pointer\n";
                    throw std::runtime_error("getEnergy return is NULL");
                }
                std::cout << "    done calling get energy\n";

                // parse the returned tuple
                double energy = PyFloat_AsDouble(returnval);

                // decrease referenece counts on Python objects
                Py_XDECREF(returnval); 

                return energy;
            }

            virtual double get_energy_gradient(Array<double> x, Array<double> grad)
            {
                assert(x.size() == grad.size());

                // create a numpy array from x
                // this wraps the data which is unsafe becase the python object might live longer than the data in x.data
                std::cout << "in get_energy_gradient\n";
                npy_intp N = (npy_intp) x.size();
                PyObject * numpyx = PyArray_SimpleNewFromData(1, &N, NPY_DOUBLE, static_cast<void*>(x.data()));
                
                // call the function getEnergy
                PyObject * name = PyString_FromString("getEnergyGradient");
                PyObject * returnval = PyObject_CallMethodObjArgs(_potential, name, numpyx, NULL);
                if (!returnval){
                    //parse error
                    throw std::runtime_error("getEnergyGradient return is NULL");
                }
                std::cout << "    done calling get energy gradient\n";

                // parse the returned tuple
                double energy;
                PyObject * numpy_grad; //the reference count for this does not need to be decreased
                if (!PyArg_ParseTuple(returnval, "dO", &energy, &numpy_grad)) 
                    throw std::runtime_error("failed to parse the tuple");
                std::cout << "    done parsing tuple\n";

                //copy the gradient into grad
                double * gdata = (double*) PyArray_DATA(numpy_grad);
                for (size_t i = 0; i < grad.size(); ++i){
                    grad[i] = gdata[i];
                }
                std::cout << "    done copying grad\n";

                // decrease referenece counts on Python objects
                Py_XDECREF(numpyx); 
                Py_XDECREF(name); 
                Py_XDECREF(returnval); 
                //Py_XDECREF(numpy_grad);

                return energy;
            }
    };
}
#endif
