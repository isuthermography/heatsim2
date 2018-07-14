// from: http://stackoverflow.com/questions/17511309/fast-string-array-cython

//#cdef char ** to_cstring_array(list_str):
//#    cdef char **ret = <char **>malloc(len(list_str) * sizeof(char *))
//#    for i in xrange(len(list_str)):
//#        ret[i] = PyString_AsString(list_str[i])
//#    return ret
#include <assert.h>

char **to_cstring_array(PyObject *list_str)
{
  size_t i;
  size_t list_len;
  
  char **ret;

  assert(PyList_Check(list_str));
  list_len = PyList_Size(list_str);
  
  ret = (char **)malloc(list_len * sizeof(char *));
  for (i=0; i < list_len; i++) {
#if PY_MAJOR_VERSION < 3
    ret[i] = PyString_AsString(PyList_GetItem(list_str,i));
#else
    ret[i] = PyBytes_AsString(PyList_GetItem(list_str,i));    
#endif
  }
  return ret;
}


static inline PyObject *CObject_Or_Capsule_New(void *Ptr)
{
  PyObject *Obj;
#ifdef Py_CAPSULE_H
  Obj=PyCapsule_New(Ptr,NULL,NULL);
#else
  Obj=PyCObject_FromVoidPtr(Ptr,NULL);
#endif
  return Obj;
}

static inline void *CObject_Or_Capsule_AsVoidPtr(PyObject *Obj)
{
#ifdef Py_CAPSULE_H
  return PyCapsule_GetPointer(Obj,NULL);
#else
  return PyCObject_AsVoidPointer(Obj);
#endif
  
}
