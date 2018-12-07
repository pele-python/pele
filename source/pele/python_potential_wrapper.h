#ifndef _PELE_POTENTIAL_FUNCTION_H
#define _PELE_POTENTIAL_FUNCTION_H
#include "array.h"
#include "base_potential.h"
#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <stdexcept>

//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

namespace pele {

/**
 * This class wraps a get_energy function and a get_energy_gradient function in
 * a class that derives from BasePotential.  This is necessary to be able to use
 * the functions in the pele c++ interface.  This is a backup method, the
 * preferred method is to define a class separately for each potential, which
 * eliminates the need to pass around the void * userdata parameter.
 */
/*
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
*/


/**
 * This class derives from the c++ BasePotential, but wraps a pure python
 * potential This is necessary to be able to use the functions in the pele
 * c++ interface.
 */
class PythonPotential : public BasePotential
{
    PyObject * _potential;

public:
    PythonPotential(PyObject * potential)
        :    _potential(potential)
    {
        Py_XINCREF(_potential);

        // import the the numpy array API.  This is commented because
        // it is now done in the cython module file.  It is possible I
        // need to define some preprocessor variables like
        // NO_IMPORT_ARRAY and PY_ARRAY_UNIQUE_SYMBOL, but the
        // documentation is a bit confusing.
        // http://docs.scipy.org/doc/numpy/reference/c-api.array.html
        //import_array();
    }

    virtual ~PythonPotential() 
    { 
        Py_XDECREF(_potential); 
    }

    /**
     * call the getEnergy method of the python potential
     */
    virtual double get_energy(Array<double> x) 
    { 
        // create a numpy array from x
        // copy the data from x because becase the python object might
        // live longer than the data in x.data
        //otherwise this alone would do PyObject * numpyx = PyArray_SimpleNewFromData(1, &N, NPY_DOUBLE, x.data());
        npy_intp N = (npy_intp) x.size();
        PyObject * numpyx = PyArray_SimpleNew(1, &N, NPY_DOUBLE);
        if (!numpyx){
            std::cerr << "created numpy object is NULL\n";
            throw std::runtime_error("created numpy object is NULL\n");
        }
        double * xdata = (double*) PyArray_DATA(numpyx);
        for (size_t i = 0; i < x.size(); ++i){
            xdata[i] = x[i];
        }


        // call the function getEnergy
#if PY_MAJOR_VERSION >= 3
        PyObject * name = PyUnicode_FromString("getEnergy");
#else
        PyObject * name = PyString_FromString("getEnergy");
#endif
        PyObject * returnval = PyObject_CallMethodObjArgs(_potential, name, numpyx, NULL);
        Py_XDECREF(name); 
        Py_XDECREF(numpyx); 
        if (!returnval){
            //parse error
            throw std::runtime_error("getEnergy returned NULL");
        }
        //std::cout << "    done calling get energy\n";

        // parse the returned tuple
        double energy = PyFloat_AsDouble(returnval);
        Py_XDECREF(returnval);
        //TODO: error checking
        if (PyErr_Occurred()){
            PyErr_Clear();
            PyErr_SetString(PyExc_TypeError, "return value of getEnergy could not be converted to float");
            throw std::runtime_error("return value of getEnergy could not be converted to float");
        }

        return energy;
    }

    /**
     * call the getEnergyGradient method of the python potential
     */
    virtual double get_energy_gradient(Array<double> x, Array<double> grad)
    {
        if (x.size() != grad.size()) {
            throw std::invalid_argument("grad.size() be the same as x.size()");
        }

        // create a numpy array from x
        // copy the data from x because because the python object might
        // live longer than the data in x.data
        npy_intp N = (npy_intp) x.size();
        PyObject * numpyx = PyArray_SimpleNew(1, &N, NPY_DOUBLE);
        if (!numpyx){
            std::cerr << "created numpy object is NULL\n";
            throw std::runtime_error("created numpy object is NULL\n");
        }
        double * numpyx_data = (double*) PyArray_DATA(numpyx);
        for (size_t i = 0; i < x.size(); ++i){
            numpyx_data[i] = x[i];
        }
        
        // call the function getEnergy
#if PY_MAJOR_VERSION >= 3
        PyObject * name = PyUnicode_FromString("getEnergyGradient");
#else
        PyObject * name = PyString_FromString("getEnergyGradient");
#endif
        PyObject * returnval = PyObject_CallMethodObjArgs(_potential, name, numpyx, NULL);
        Py_XDECREF(numpyx); 
        Py_XDECREF(name); 
        if (!returnval){
            //parse error
            throw std::runtime_error("getEnergyGradient return is NULL");
        }

        // parse the returned tuple into a doulbe and a numpy array
        double energy;
        PyObject * npgrad_returned; //the reference count for this does not need to be decreased
        if (!PyArg_ParseTuple(returnval, "dO", &energy, &npgrad_returned)){
            Py_XDECREF(returnval);
            throw std::runtime_error("failed to parse the tuple");
        }

        // convert the returned gradient into an array which I know I
        // can safely use as a double array.
        // note: NPY_CARRAY is for numpy version 1.6, for later version use NPY_ARRAY_CARRAY
        PyObject * npgrad_safe = PyArray_FromAny(npgrad_returned,
                PyArray_DescrFromType(NPY_DOUBLE), 1, 1, NPY_CARRAY,
                NULL);
        if (!npgrad_safe){
            Py_XDECREF(returnval);
            throw std::runtime_error("gradient returned by getEnergyGradient could not be converted to numpy double array");
        }
        // check the size of the gradient array
        if (static_cast<size_t>(PyArray_Size(npgrad_safe)) != grad.size()){
            PyErr_SetString(PyExc_IndexError, "gradient returned by getEnergyGradient has wrong size.");
            Py_XDECREF(returnval);
            Py_XDECREF(npgrad_safe);
            throw std::runtime_error("return value of getEnergy could not be converted to float");
        }

        //copy the gradient into grad
        double * gdata = (double*) PyArray_DATA(npgrad_safe);
        for (size_t i = 0; i < grad.size(); ++i){
            grad[i] = gdata[i];
        }
        //std::cout << "    done copying grad\n";

        // decrease referenece counts on Python objects
        Py_XDECREF(returnval);
        Py_XDECREF(npgrad_safe);

        return energy;
    }
};
}
#endif
