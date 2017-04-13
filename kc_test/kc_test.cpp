// kc_test.cpp : Defines the initialization routines for the DLL.
//

#include "stdafx.h"
#include <afxwin.h>
#include <afxdllx.h>

static AFX_EXTENSION_MODULE kc_test_DLL = { NULL, NULL };

extern "C" int APIENTRY
DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID lpReserved)
    {
    UNREFERENCED_PARAMETER(lpReserved);

    if (dwReason == DLL_PROCESS_ATTACH)
        {
        TRACE0("kc_test.DLL Initializing!\n");
        
        if (!AfxInitExtensionModule(kc_test_DLL, hInstance))
            return(0);
        }
    else if (dwReason == DLL_PROCESS_DETACH)
        {
        TRACE0("kc_test.DLL Terminating!\n");
        AfxTermExtensionModule(kc_test_DLL);
        }

    return(1);
    }
