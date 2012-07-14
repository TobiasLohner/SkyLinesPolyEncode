#include <Python.h>
#include "structmember.h"
#include <list>

#include "SkyLinesPolyEncoder.h"

using namespace std;

typedef struct {
    PyObject_HEAD
    SkyLinesPolyEncoder* spe;
} SkyLinesPolyEncoderPy;

static int
SkyLinesPolyEncoderPy_init(SkyLinesPolyEncoderPy *self, PyObject *args, PyObject *kwargs)
{
    // FIXME: how to specify method signature?
    // defaults
    int numLevels=18;
    int zoomFactor=2;
    double threshold=0.00001;
    int forceEndpoints=1;
    
    static char *kwlist[] = {"num_levels", "zoom_factor", "threshold", "force_endpoints", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|iidi", kwlist, 
                                      &numLevels, &zoomFactor, &threshold, &forceEndpoints))
        return NULL; 

    self->spe = new SkyLinesPolyEncoder(numLevels, zoomFactor, threshold, forceEndpoints);
    return 0;
}

static void
SkyLinesPolyEncoderPy_dealloc(SkyLinesPolyEncoderPy *self) {
    delete self->spe;
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject*
SkyLinesPolyEncoderPy_classify(SkyLinesPolyEncoderPy *self, PyObject *args, PyObject *kwargs) {
    // FIXME: how to specify method signature?
    PyObject *p_points;
    char *type = "\0";
    bool remove = true;
   
    static char *kwlist[] = {"p_points", "remove", "type", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|bs", kwlist,
                                      &p_points, &remove, &type))
        return NULL;

    // Check we actually have a sequence of two-tuples and populate a
    // std::vector with the values. It'd be nice to pass the values through
    // without iterating needlessly...
    if (!(p_points = PySequence_Fast(p_points, "expected sequence type"))) {
        return NULL;
    }
    Py_ssize_t nPoints = PySequence_Fast_GET_SIZE(p_points);
    vector<vector<double> > n_points;
    for (Py_ssize_t i=0; i<nPoints; i++) {
        // check the item is a tuple
        PyObject *p_vertices = PySequence_Fast(PySequence_Fast_GET_ITEM(p_points, i), "expected sequence type");
        if (!p_vertices) {
            return NULL;
        }
        // and its length is at least 2
        if (PySequence_Fast_GET_SIZE(p_vertices) < 2) {
            PyErr_SetString(PyExc_TypeError, "expected (n>=2)-tuple with numbers ((x0, y0), (x1, y1), ...)");
            return NULL;
        } 
        // get the first two values from the item and check they're numbers
        vector<double> n_c;
        for (int j=0; j<PySequence_Fast_GET_SIZE(p_vertices); j++) {
            PyObject *p_c = PySequence_Fast_GET_ITEM(p_vertices, j);
            if (!PyNumber_Check(p_c)) {
                PyErr_SetString(PyExc_TypeError, "expected (n>=2)-tuple with numbers ((x0, y0), (x1, y1), ...)");
                return NULL;
            }
            // get as double
            PyObject *pf_c = PyNumber_Float(p_c);
            n_c.push_back(PyFloat_AsDouble(pf_c));
            Py_DECREF(pf_c);
        }
        // add to our vector
        n_points.push_back(n_c);
        Py_DECREF(p_vertices);
    }
    Py_DECREF(p_points);
    
    // do our encoding
    vector<int> n_levels;
    // might as well allow some other threads to run...
    Py_BEGIN_ALLOW_THREADS
    n_levels = self->spe->dpEncode(n_points, type);
    Py_END_ALLOW_THREADS
    
    // build a dictionary of the results
    // {'points':[...], 'levels':[...], 'numLevels':#, 'zoomFactor':#}
    PyObject *p_result = PyDict_New();
    PyObject *ep, *el, *zf, *nl;
    
    PyObject *points = PyList_New(0);
    PyObject *levels = PyList_New(0);

    for (size_t i=0; i < n_points.size(); i++) {
      if (remove && n_levels[i] == -1) continue;

      ep = PyTuple_New(n_points.at(i).size());
      int j = 0;
      for (vector<double>::iterator it = n_points.at(i).begin(); it != n_points.at(i).end(); it++) {
        PyTuple_SetItem(ep, j, PyFloat_FromDouble( *it ));
        j++;
      }

      if (PyList_Append(points, ep))
        return NULL;

    }
    Py_DECREF(ep);

    if (PyDict_SetItemString(p_result, "points", points))
      return NULL;
    Py_DECREF(points);


    for (size_t i=0; i < n_levels.size(); i++) {
      if (remove && n_levels[i] == -1) continue;
      
      el = PyInt_FromLong(n_levels[i]);
      if (PyList_Append(levels, el))
        return NULL;
    }
    Py_DECREF(el);

    if (PyDict_SetItemString(p_result, "levels", levels))
      return NULL;
    Py_DECREF(levels);

    zf = PyInt_FromLong(self->spe->getZoomFactor());
    if (PyDict_SetItemString(p_result, "zoomFactor", zf))
        return NULL;
    Py_DECREF(zf);
    
    nl = PyInt_FromLong(self->spe->getNumLevels());
    if (PyDict_SetItemString(p_result, "numLevels", nl))
        return NULL;
    Py_DECREF(nl);
    
    return p_result;
}

static PyObject*
SkyLinesPolyEncoderPy_encodeList(SkyLinesPolyEncoderPy *self, PyObject *args) {
    // FIXME: how to specify method signature?
    PyObject *p_list;

    if (!PyArg_ParseTuple(args, "O", &p_list))
        return NULL;

    // Check we actually have a sequence of two-tuples and populate a
    // std::vector with the values. It'd be nice to pass the values through
    // without iterating needlessly...
    if (!(p_list = PySequence_Fast(p_list, "expected sequence type"))) {
        return NULL;
    }

    Py_ssize_t nList = PySequence_Fast_GET_SIZE(p_list);

    list<int> n_list;
    for (Py_ssize_t i=0; i<nList; i++) {
        PyObject *p_c = PySequence_Fast_GET_ITEM(p_list, i);
        if (!PyNumber_Check(p_c)) {
          PyErr_SetString(PyExc_TypeError, "expected list with numbers");
          return NULL;
        }

        PyObject *pn = PyNumber_Int(p_c);
        n_list.push_back(PyInt_AsLong(pn));
        Py_DECREF(pn);
    }
    Py_DECREF(p_list);

    // do our encoding
    string n_result;
    // might as well allow some other threads to run...
    Py_BEGIN_ALLOW_THREADS
    n_result = self->spe->encodeList(n_list);
    Py_END_ALLOW_THREADS

    // return the result
    PyObject *el;

    el = PyString_FromString(n_result.c_str());

    return el;
}

static PyObject*
SkyLinesPolyEncoderPy_encode(SkyLinesPolyEncoderPy *self, PyObject *args) {
    // FIXME: how to specify method signature?
    PyObject *p_points;
    PyObject *p_levels;
    if (!PyArg_ParseTuple(args, "OO", &p_points, &p_levels))
        return NULL;

    // Check we actually have a sequence of two-tuples and populate a
    // std::vector with the values. It'd be nice to pass the values through
    // without iterating needlessly...
    if (!(p_points = PySequence_Fast(p_points, "expected sequence type"))) {
        return NULL;
    }

    if (!(p_levels = PySequence_Fast(p_levels, "expected sequence type"))) {
        return NULL;
    }

    Py_ssize_t nPoints = PySequence_Fast_GET_SIZE(p_points);
    Py_ssize_t nLevels = PySequence_Fast_GET_SIZE(p_levels);

    if (nPoints != nLevels) {
      PyErr_SetString(PyExc_ValueError, "number of points not equal number of levels");
      return NULL;
    }

    vector<pair<double,double> > n_points;
    vector<int> n_levels;
    for (Py_ssize_t i=0; i<nPoints; i++) {
        // check the item is a tuple
        PyObject *p_vertices = PySequence_Fast(PySequence_Fast_GET_ITEM(p_points, i), "expected sequence type");
        if (!p_vertices) {
            return NULL;
        }
        // and its length is at least 2 (ignore 3D coordinates if present)
        if (PySequence_Fast_GET_SIZE(p_vertices) < 2) {
            PyErr_SetString(PyExc_TypeError, "expected two-tuple with numbers ((x0, y0), (x1, y1), ...)");
            return NULL;
        }
        // get the first two values from the item and check they're numbers
        double n_c[2];
        for (int j=0; j<2; j++) {
            PyObject *p_c = PySequence_Fast_GET_ITEM(p_vertices, j);
            if (!PyNumber_Check(p_c)) {
                PyErr_SetString(PyExc_TypeError, "expected two-tuple with numbers ((x0, y0), (x1, y1), ...)");
                return NULL;
            }
            // get as double
            PyObject *pf_c = PyNumber_Float(p_c);
            n_c[j] = PyFloat_AsDouble(pf_c);
            Py_DECREF(pf_c);
        }
        // add to our vector
        n_points.push_back(pair<double,double>(n_c[0], n_c[1]));
        Py_DECREF(p_vertices);

        PyObject *p_c = PySequence_Fast_GET_ITEM(p_levels, i);
        if (!PyNumber_Check(p_c)) {
          PyErr_SetString(PyExc_TypeError, "expected list with numbers");
          return NULL;
        }

        PyObject *pn = PyNumber_Int(p_c);
        n_levels.push_back(PyInt_AsLong(pn));
        Py_DECREF(pn);
    }
    Py_DECREF(p_points);
    Py_DECREF(p_levels);

    // do our encoding
    auto_ptr<pair<string, string> > n_result;
    // might as well allow some other threads to run...
    Py_BEGIN_ALLOW_THREADS
    n_result = self->spe->encode(n_points, n_levels);
    Py_END_ALLOW_THREADS

    // build a dictionary of the results
    // {'points':'...', 'levels':'...'}
    PyObject *p_result = PyDict_New();
    PyObject *ep, *el;

    ep = PyString_FromString(n_result->first.c_str());
    if (PyDict_SetItemString(p_result, "points", ep))
        return NULL;
    Py_DECREF(ep);

    el = PyString_FromString(n_result->second.c_str());
    if (PyDict_SetItemString(p_result, "levels", el))
        return NULL;
    Py_DECREF(el);

    return p_result;
}


static PyMemberDef SkyLinesPolyEncoderPy_members[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef SkyLinesPolyEncoderPy_methods[] = {
    {"classify", (PyCFunction)SkyLinesPolyEncoderPy_classify, METH_VARARGS|METH_KEYWORDS, "Classify a sequence of points"},
    {"encode", (PyCFunction)SkyLinesPolyEncoderPy_encode, METH_VARARGS, "Encode a sequence of classified (lon,lat)-points"},
    {"encodeList", (PyCFunction)SkyLinesPolyEncoderPy_encodeList, METH_VARARGS, "Encode a sequence of numbers"},
    {NULL}  /* Sentinel */
};

static PyTypeObject SkyLinesPolyEncoderPyType = {
    PyObject_HEAD_INIT(NULL)
    0,                                    /*ob_size*/
    "skylinespolyencode.SkyLinesPolyEncoder", /*tp_name*/
    sizeof(SkyLinesPolyEncoderPy),        /*tp_basicsize*/
    0,                                    /*tp_itemsize*/
    (destructor)SkyLinesPolyEncoderPy_dealloc,   /*tp_dealloc*/
    0,                                    /*tp_print*/
    0,                                    /*tp_getattr*/
    0,                                    /*tp_setattr*/
    0,                                    /*tp_compare*/
    0,                                    /*tp_repr*/
    0,                                    /*tp_as_number*/
    0,                                    /*tp_as_sequence*/
    0,                                    /*tp_as_mapping*/
    0,                                    /*tp_hash */
    0,                                    /*tp_call*/
    0,                                    /*tp_str*/
    0,                                    /*tp_getattro*/
    0,                                    /*tp_setattro*/
    0,                                    /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,       /*tp_flags*/
    "SkyLines Polyline Encoder",       /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    SkyLinesPolyEncoderPy_methods,               /* tp_methods */
    SkyLinesPolyEncoderPy_members,               /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    (initproc)SkyLinesPolyEncoderPy_init,        /* tp_init */
    0,                                    /* tp_alloc */
    PyType_GenericNew,                    /* tp_new */
};

static PyMethodDef skylinespolyencode_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initskylinespolyencode(void) 
{
    PyObject* m;

    if (PyType_Ready(&SkyLinesPolyEncoderPyType) < 0)
        return;

    m = Py_InitModule3("skylinespolyencode", skylinespolyencode_methods,
                       "SkyLines Polyline encoding (C++ extension)");
                       
    if (m == NULL)
        return;

    Py_INCREF(&SkyLinesPolyEncoderPyType);
    PyModule_AddObject(m, "SkyLinesPolyEncoder", (PyObject *)&SkyLinesPolyEncoderPyType);
}

