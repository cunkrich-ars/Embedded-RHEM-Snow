// snow.cpp : Defines the functions for the static library.

#include "stdlib.h"
#include "pch.h"
#include "framework.h"
#include <C:\Users\Carl\AppData\Local\Programs\Python\Python311\include\Python.h>

static PyObject *pModule;
static PyObject *pFunc_get_next_event;
static PyObject *pFunc_get_npoints;
static PyObject *pFunc_get_times;
static PyObject *pFunc_get_depths;
static PyObject *pFunc_get_sat;
static PyObject *pFunc_get_ice;
static PyObject *pList;
static PyObject *pItem;
static PyObject *pValue;

static long npoints;


extern "C" long snow_run(char *ids, char *forcing_files, char *soils, double slopes, double aspects) // use with Fortran
//long snow_run(char* ids, char* forcing_files, char* soils, double slopes, double aspects)
{
    long flag;
    PyConfig config;
    PyObject *pName;
    PyObject *pFunc_run;
    PyObject *pArgs;
    PyStatus status;

    FILE* fp;
 //   FILE* fp2;
    char path_str[10][200];
    wchar_t w_path_str[10][200];
    size_t nchars;
//    wchar_t wc;
//    wchar_t *wstr;
    errno_t err;
    int ierr, j;
//    int i;

    Py_Initialize();

    ierr = Py_IsInitialized();
    if (ierr == 0) return 0;

    PyConfig_InitPythonConfig(&config);
    config.module_search_paths_set = 1;

    // append local dir to path

    PyWideStringList_Append(&config.module_search_paths, L".");

    // Read file with module search paths

    err = fopen_s(&fp, "modpaths.txt", "r");
 //   err = fopen_s(&fp2, "paths.txt", "w");
 //   ierr = fwide(fp2, 1);
    j = 0;
//    while (j < 5) {
    while (fgets(path_str[j], 200, fp) != NULL) {
        path_str[j][strcspn(path_str[j], "\n")] = 0;
        err = mbstowcs_s(&nchars, w_path_str[j], 200, path_str[j], 199);
        status = PyWideStringList_Append(&config.module_search_paths, w_path_str[j]);
        if (PyStatus_Exception(status)) return 1;
//        for (i = 0; w_path_str[j][i] != L'\0'; i++) {
//            wc = putwchar(w_path_str[j][i]);
//            wc = putwc(w_path_str[j][i], fp2);
//        }
        j++;
    }
    fclose(fp);
//    fclose(fp2);

//    PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\Carl\\AppData\\Local\\Programs\\Python\\Python311\\Lib");
//    PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\Carl\\AppData\\Local\\Programs\\Python\\Python311\\python311.zip");
//    PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\Carl\\AppData\\Local\\Programs\\Python\\Python311\\DLLs");
//    PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\Carl\\AppData\\Local\\Programs\\Python\\Python311");
//    PyWideStringList_Append(&config.module_search_paths, L"C:\\Users\\Carl\\AppData\\Local\\Programs\\Python\\Python311\\Lib\\site-packages");

    status = Py_InitializeFromConfig(&config);

    if (PyStatus_Exception(status)) return 2;
 
    pName = PyUnicode_DecodeFSDefault("snow");
    pModule = PyImport_Import(pName);
 
    if (pModule == NULL) return 3;
 
//    ierr = puts("snow module imported\n");

    pArgs = PyTuple_New(5); // argument list passed to python

// convert each argument to a list object and add to pArgs

    pList = PyList_New(1);
    pItem = PyUnicode_DecodeFSDefault(ids);
    PyList_SetItem(pList, 0, pItem);
    PyTuple_SetItem(pArgs, 0, pList);

    pList = PyList_New(1);
    pItem = PyUnicode_DecodeFSDefault(forcing_files);
    PyList_SetItem(pList, 0, pItem);
    PyTuple_SetItem(pArgs, 1, pList);

    pList = PyList_New(1);
    pItem = PyUnicode_DecodeFSDefault(soils);
    PyList_SetItem(pList, 0, pItem);
    PyTuple_SetItem(pArgs, 2, pList);

    pList = PyList_New(1);
    pItem = PyFloat_FromDouble(slopes);
    PyList_SetItem(pList, 0, pItem);
    PyTuple_SetItem(pArgs, 3, pList);

    pList = PyList_New(1);
    pItem = PyFloat_FromDouble(aspects);
    PyList_SetItem(pList, 0, pItem);
    PyTuple_SetItem(pArgs, 4, pList);

    pFunc_run = PyObject_GetAttrString(pModule, "run");

    if (pFunc_run == NULL) return 5;

//    ierr = puts("calling run\n");

    pValue = PyObject_CallObject(pFunc_run, pArgs);

    if (pValue != NULL) {
        flag = PyLong_AsLong(pValue);
    }
    else {
        return 6;
    }
 
    pFunc_get_npoints = PyObject_GetAttrString(pModule, "get_npoints");
    pFunc_get_next_event = PyObject_GetAttrString(pModule, "get_next_event");
    pFunc_get_times = PyObject_GetAttrString(pModule, "get_times");
    pFunc_get_depths = PyObject_GetAttrString(pModule, "get_depths");
    pFunc_get_sat = PyObject_GetAttrString(pModule, "get_sat");
    pFunc_get_ice = PyObject_GetAttrString(pModule, "get_ice");

    Py_DECREF(pName);
    Py_DECREF(pArgs);
    Py_DECREF(pFunc_run);

    return flag;
}

extern "C" long snow_next_event(long *date) // use with Fortran
//long snow_next_event(long *date)
{
    long i;

    pList = PyObject_CallObject(pFunc_get_next_event, NULL);

    for (i = 0; i < 3; i++)
    {
        pItem = PyList_GetItem(pList, i);
        date[i] = PyLong_AsLong(pItem);
    }

    return 0;
}

extern "C" long snow_npoints(void) // use with Fortran
//long snow_npoints(void)
{
    pValue = PyObject_CallObject(pFunc_get_npoints, NULL);
    npoints = PyLong_AsLong(pValue);

    return npoints;
}

extern "C" long snow_times(double *t) // use with Fortran
//long snow_times(double *t)
{
    long i;

    pList = PyObject_CallObject(pFunc_get_times, NULL);

    for (i = 0; i < npoints; i++)
    {
        pItem = PyList_GetItem(pList, i);
        t[i] = PyFloat_AsDouble(pItem);
    }

    return 0;
}

extern "C" long snow_depths(double *d) // use with Fortran
//long snow_depths(double *d)
{
    long i;

     pList = PyObject_CallObject(pFunc_get_depths, NULL);

    for (i = 0; i < npoints; i++)
    {
        pItem = PyList_GetItem(pList, i);
        d[i] = PyFloat_AsDouble(pItem);
    }

    return 0;
}

extern "C" double snow_sat(void) // use with Fortran
//double snow_sat(void)
{
    double sat;

    pValue = PyObject_CallObject(pFunc_get_sat, NULL);
    sat = PyFloat_AsDouble(pValue);
 
    return sat;
}

extern "C" double snow_ice(void) // use with Fortran
//double snow_ice(void)
{
    double ice;

    pValue = PyObject_CallObject(pFunc_get_ice, NULL);
    ice = PyFloat_AsDouble(pValue);

    return ice;
}

extern "C" long snow_finale(void) // use with Fortran
//long snow_finale(void)
{
    Py_DECREF(pModule);
    Py_DECREF(pValue);
    Py_DECREF(pList);
    Py_DECREF(pItem);
    Py_DECREF(pFunc_get_next_event);
    Py_DECREF(pFunc_get_npoints);
    Py_DECREF(pFunc_get_times);
    Py_DECREF(pFunc_get_depths);
    Py_DECREF(pFunc_get_sat);
    Py_DECREF(pFunc_get_ice);

    if (Py_FinalizeEx() < 0) {
        return 1;
    }
    return 0;
}
