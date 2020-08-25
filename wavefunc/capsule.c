#define PY_SSIZE_T_CLEAN
#include <Python.h>

void capsule_cleanup(PyObject *capsule) {
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    // I'm going to assume your memory needs to be freed with free().
    // If it needs different cleanup, perform whatever that cleanup is
    // instead of calling free().
    free(memory);
}