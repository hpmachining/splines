#include "CDEHead.h"

// User Include Section
#include "Funcs.h"


// Data Section
static const CDE_data CDE_FUN_TestSplineLibrary[] =
    {
    _T("TestSplineLibrary"),        // name
    NULL,
    1,                              // # of parameters
        _T("TestSplineLibrary"),    // function name
        _T(""),                     // array size
        CDE_INT,                    // return value type
    };


static const CDE_data *Data[] =
    {
    CDE_FUN_TestSplineLibrary,
    };



// Function Section
static const CDE_code CDE_TestSplineLibrary = 
    {
    _T("TestSplineLibrary"),        // Function name
    NULL,                           // (Not hidden)
    NULL,                           // (No alias)
    CDE_ANY,                        // Flags
    CDE_FUN_TestSplineLibrary,      // Parameter list
    (int (*)())TestSplineLibrary,   // Function Pointer
    };


static const CDE_code *Functions[] = {
    &CDE_TestSplineLibrary,
};


// App Section

// Popup Section

// Import Section

// Export Section

// Unassociated Section

// Classic Menu Section

// Main Section
extern "C" __declspec(dllexport) const CDE_module CDEMODS2_Funcs = 
    {
    1,                              // Number of functions
    Functions,                      // Function decl. section
    1,                              // Number of data types
    Data,                           // Data declaration section
    0,                              // Number of Applications
    NULL,                           // No Applications
    0,                              // Number of Popups
    NULL,                           // No Popups
    0,                              // Number of Imports
    NULL,                           // No Imports
    0,                              // Number of Exports
    NULL,                           // No Exports
    0,                              // Number of Unassociated Functions
    NULL,                           // No Unassociated Functions
    0,                              // Number of Classic Menus
    NULL,                           // No Classic Menus
    NULL,                           // No Help File Name
    NULL,                           // No Resource DLL Name
    NULL,                           // No Product Name
    };