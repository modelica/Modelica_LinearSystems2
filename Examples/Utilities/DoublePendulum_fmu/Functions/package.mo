within Modelica_LinearSystems2.Examples.Utilities.DoublePendulum_fmu;
package Functions "DoublePendulum FMI functions"
  extends FMI.Interfaces.PartialFmiFunctions(
nx = 6,
nu = 1,
ny = 6,
  id_u={352321536},
  id_y={335544325,335544324,335544320,335544321,335544323,335544322});


  redeclare class FmiInstance
    extends ExternalObject;
    function constructor "Initialize FMI model"
      extends Modelica.Icons.Function;
      input String instanceName;
      input Boolean loggingOn;
      output FmiInstance fmi;
    external"C" fmi = DoublePendulum_fmiInstantiateModel2(instanceName,
        loggingOn)
    annotation(Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_Instantiate_C
#define DoublePendulum_Instantiate_C 1
#include <stdlib.h>
void DoublePendulumLogger(fmiComponent c, fmiString instanceName, fmiStatus status,
                           fmiString category, fmiString message, ...)  { 
  char msg[4000];
  char buf[4000];
  int len;
  va_list ap;
  va_start(ap,message);
#if defined(_MSC_VER) && _MSC_VER>=1200
  len = _snprintf(msg, sizeof(msg)/sizeof(*msg), \"%s: %s\", instanceName, message);
  if (len < 0) goto fail;
  len = _vsnprintf(buf, sizeof(buf)/sizeof(*buf) - 2, msg, ap);
  if (len < 0) goto fail;
#else
  len = vsprintf(buf, message, ap);
#endif
  buf[len] = 0;
  va_end(ap);
  switch (status)  { 
    case fmiFatal:
      ModelicaError(buf);
      break;
    default:
      ModelicaMessage(buf);
      break;
   } 
  return;
fail:
  ModelicaMessage(\"Logger failed, message too long?\");
 } 
void * DoublePendulum_fmiInstantiateModel2(const char*instanceName, fmiBoolean loggingOn)  { 
  static fmiCallbackFunctions funcs =  { &DoublePendulumLogger, &calloc, &free } ;
  struct dyDoublePendulum_Extended* res;

  res = calloc(1, sizeof(struct dyDoublePendulum_Extended));
  if (res!=0)  { 
    res->m=DoublePendulum_fmiInstantiateModel(instanceName, \" { ecb6e32d-81ac-4212-92ab-ac188d5f98ce } \", funcs, loggingOn);
    if (0==res->m)  { free(res);res=0;ModelicaError(\"InstantiateModel failed\"); } 
    else  { res->dyTriggered=0;res->dyTime=res->dyLastTime=-1e37; } 
   } 
  return res;
 } 
#endif",
     Library="DoublePendulum");
    end constructor;

    function destructor "Release storage of FMI model"
      extends Modelica.Icons.Function;
      input FmiInstance fmi;
    external"C" DoublePendulum_fmiFreeModelInstance2(fmi);
      annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_Free_C
#define DoublePendulum_Free_C 1
#include <stdlib.h>
void DoublePendulum_fmiFreeModelInstance2(void*m)  { 
  struct dyDoublePendulum_Extended*a=m;
  if (a)  { 
    DoublePendulum_fmiFreeModelInstance(a->m);
    free(a);
   } 
 } 
#endif", Library="DoublePendulum");
    end destructor;
  end FmiInstance;


  redeclare function extends fmiSetTime
  external"C" DoublePendulum_fmiSetTime2(fmi, ti);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_SetTime_C
#define DoublePendulum_SetTime_C 1
#include <stdlib.h>
void DoublePendulum_fmiSetTime2(void*m, double ti)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
    a->dyTime=ti;
    status=DoublePendulum_fmiSetTime(a->m, ti);
   } 
  if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"SetTime failed\");
 } 
#endif", Library="DoublePendulum");
  end fmiSetTime;


  redeclare function extends fmiSetContinuousStates
  external"C" DoublePendulum_fmiSetContinuousStates2(
        fmi,
        x,
        size(x, 1));
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_SetContinuousStates_C
#define DoublePendulum_SetContinuousStates_C 1
#include <stdlib.h>
void DoublePendulum_fmiSetContinuousStates2(void*m, const double*x, size_t nx)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiSetContinuousStates(a->m, x, nx);
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"SetContinuousStates failed\");
 } 
#endif", Library="DoublePendulum");
  end fmiSetContinuousStates;


  redeclare function extends fmiGetContinuousStates
  external"C" DoublePendulum_fmiGetContinuousStates2(
        fmi,
        x,
        nx);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_GetContinuousStates_C
#define DoublePendulum_GetContinuousStates_C 1
#include <stdlib.h>
void DoublePendulum_fmiGetContinuousStates2(void*m, double*x, int nx)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiGetContinuousStates(a->m, x, nx);
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"GetContinuousStates failed\");
 } 
#endif", Library="DoublePendulum");
  end fmiGetContinuousStates;


  redeclare function extends fmiCompletedStep
  external"C" crossing=  DoublePendulum_fmiCompletedStep2(fmi);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_CompletedStep_C
#define DoublePendulum_CompletedStep_C 1
#include <stdlib.h>
double DoublePendulum_fmiCompletedStep2(void*m)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     if (a->dyTime>a->dyLastTime)  { 
       fmiBoolean b=0;
       status=DoublePendulum_fmiCompletedIntegratorStep(a->m, &b);
       a->dyLastTime=a->dyTime;
       if (b) a->dyTriggered=1;
      }  else status=fmiOK;
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"CompletedIntegratorStep failed\");
   return a->dyTriggered && a->dyTime>=a->dyLastTime;
 } 
#endif", Library="DoublePendulum");
  end fmiCompletedStep;


  redeclare function extends fmiEventUpdate
  external"C" stateReset=  DoublePendulum_fmiEventUpdate2(fmi, tnext)
  annotation(Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_EventUpdate_C
#define DoublePendulum_EventUpdate_C 1
#include <stdlib.h>
int DoublePendulum_fmiEventUpdate2(void*m, double*tnext) 
 { 
  struct dyDoublePendulum_Extended*a=m;
    fmiEventInfo ev;
    fmiStatus status=fmiFatal;
    ev.nextEventTime=1e37;
    if (a)  { 
    status=DoublePendulum_fmiEventUpdate(a->m, 0, &ev);
     a->dyTriggered=0;
     a->dyLastTime=a->dyTime;
    } 
   if (ev.terminateSimulation) terminate(\"Terminate signaled by FMU\");
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"EventUpdate failed\");
  *tnext=ev.nextEventTime;
  return ev.stateValuesChanged;
 } 
#endif",
     Library="DoublePendulum");
  end fmiEventUpdate;


  redeclare function extends fmiInitialize
  external"C" tnext=  DoublePendulum_fmiInitialize2(fmi);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_Initialize_C
#define DoublePendulum_Initialize_C 1
#include <stdlib.h>
double DoublePendulum_fmiInitialize2(void*m)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  fmiEventInfo ev;
  ev.nextEventTime=1e37;
  if (a)  { 
     status=DoublePendulum_fmiInitialize(a->m, fmiFalse, 1e-4, &ev);
     a->dyTriggered=0;
     a->dyLastTime=a->dyTime;
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"Initialize failed\");
   return ev.nextEventTime;
 } 
#endif", Library="DoublePendulum");
  end fmiInitialize;


  redeclare function extends fmiSetReal
  external"C" DoublePendulum_fmiSetReal2(
        fmi,
        refs,
        size(refs, 1),
        vals);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_SetReal_C
#define DoublePendulum_SetReal_C 1
#include <stdlib.h>
void DoublePendulum_fmiSetReal2(void*m, const int*refs, size_t nrefs, const double*vals)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiSetReal(a->m, refs, nrefs, vals);
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"SetReal failed\");
 } 
#endif", Library="DoublePendulum");
  end fmiSetReal;


   redeclare function extends fmiGetReal
external"C"   DoublePendulum_fmiGetReal2(
        fmi,
        refs,
        size(refs, 1),
        vals);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_GetReal_C
#define DoublePendulum_GetReal_C 1
#include <stdlib.h>
void DoublePendulum_fmiGetReal2(void*m, const int*refs, size_t nrefs, double*vals)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiGetReal(a->m, refs, nrefs, vals);
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"GetReal failed\");
 } 
#endif", Library="DoublePendulum");
   end fmiGetReal;


 redeclare function extends fmiGetInteger
  external"C" DoublePendulum_fmiGetInteger2(
        fmi,
        refs,
        size(refs, 1),
        vals);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_GetInteger_C
#define DoublePendulum_GetInteger_C 1
#include <stdlib.h>
void DoublePendulum_fmiGetInteger2(void*m, const int*refs, size_t nrefs, int*vals)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiGetInteger(a->m, refs, nrefs, vals);
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"GetInteger failed\");
 } 
#endif", Library="DoublePendulum");
 end fmiGetInteger;

//  redeclare function extends fmiSetInteger
//   external"C" DoublePendulum_fmiSetInteger2(
//         fmi,
//         refs,
//         size(refs, 1),
//         vals);
//     annotation (Header="
// #ifndef DoublePendulumFMI_C
// #define DoublePendulumFMI_C 1
// #undef MODEL_IDENTIFIER
// #define MODEL_IDENTIFIER DoublePendulum
// #include \"FMI/fmiModelFunctions.h\"
// struct dyDoublePendulum_Extended
//  { 
//   fmiComponent m;
//   double dyTime;
//   double dyLastTime;
//   int dyTriggered;
//  } ;
// #endif
// #ifndef DoublePendulum_SetInteger_C
// #define DoublePendulum_SetInteger_C 1
// #include <stdlib.h>
// void DoublePendulum_fmiSetInteger2(void*m, const int*refs, size_t nrefs, int*vals)  { 
//   struct dyDoublePendulum_Extended*a=m;
//   fmiStatus status=fmiFatal;
//   if (a)  { 
//      status=DoublePendulum_fmiSetInteger(a->m, refs, nrefs, vals);
//     } 
//    if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"SetInteger failed\");
//  } 
// #endif", Library="DoublePendulum");
//   end fmiSetInteger;


  redeclare function extends fmiGetBoolean
   external"C" DoublePendulum_fmiGetBoolean2(
        fmi,
        refs,
        size(refs, 1),
        vals);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_GetBoolean_C
#define DoublePendulum_GetBoolean_C 1
void DoublePendulum_fmiGetBoolean2(void*m, const int* refs, size_t nr, int* vals)  { 
  int i;
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiGetBoolean(a->m, refs, nr, (fmiBoolean*)(vals));
   } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"GetBoolean failed\");
   for(i=nr-1;i>=0;i--) vals[i]=((fmiBoolean*)(vals))[i];
  } 
#endif", Library="DoublePendulum");
  end fmiGetBoolean;


    redeclare function extends fmiSetBoolean
external"C"   DoublePendulum_fmiSetBoolean2(
        fmi,
        refs,
        size(refs, 1),
        vals,
        dummy);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_SetBoolean_C
#define DoublePendulum_SetBoolean_C 1
void DoublePendulum_fmiSetBoolean2(void*m, const int* refs, size_t nr, const int* vals,int*dummy)  { 
  int i;
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  for(i=0;i<nr;++i) ((fmiBoolean*)(dummy))[i]=vals[i];
  if (a)  { 
     status=DoublePendulum_fmiSetBoolean(a->m, refs, nr, (fmiBoolean*)(dummy));
   } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"SetBoolean failed\");
  } 
#endif", Library="DoublePendulum");
    end fmiSetBoolean;


  redeclare function extends fmiGetDerivatives
  external"C" DoublePendulum_fmiGetDerivatives2(
        fmi,
        dx,
        nx);
    annotation (Header="
#ifndef DoublePendulumFMI_C
#define DoublePendulumFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER DoublePendulum
#include \"FMI/fmiModelFunctions.h\"
struct dyDoublePendulum_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef DoublePendulum_GetDerivatives_C
#define DoublePendulum_GetDerivatives_C 1
#include <stdlib.h>
void DoublePendulum_fmiGetDerivatives2(void*m,double*dx,int nx)  { 
  struct dyDoublePendulum_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=DoublePendulum_fmiGetDerivatives(a->m, dx, nx);
    } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"GetDerivatives failed\");
 } 
#endif", Library="DoublePendulum");
  end fmiGetDerivatives;


  redeclare function extends fmiSetInteger
  external"C" mESCfixedT_fmiGetBoolean2(
        fmi,
        refs,
        size(refs, 1),
        vals);
    annotation (Header="
#ifndef mESCfixedTFMI_C
#define mESCfixedTFMI_C 1
#undef MODEL_IDENTIFIER 
#define MODEL_IDENTIFIER mESCfixedT
#include \"FMI/fmiModelFunctions.h\"
struct dymESCfixedT_Extended 
 { 
  fmiComponent m;
  double dyTime;
  double dyLastTime;
  int dyTriggered;
 } ;
#endif
#ifndef mESCfixedT_GetBoolean_C
#define mESCfixedT_GetBoolean_C 1
void mESCfixedT_fmiGetBoolean2(void*m, const int* refs, size_t nr, int* vals)  { 
  int i;
  struct dymESCfixedT_Extended*a=m;
  fmiStatus status=fmiFatal;
  if (a)  { 
     status=mESCfixedT_fmiGetBoolean(a->m, refs, nr, (fmiBoolean*)(vals));
   } 
   if (status!=fmiOK && status!=fmiWarning) ModelicaError(\"GetBoolean failed\");
   for(i=nr-1;i>=0;i--) vals[i]=((fmiBoolean*)(vals))[i];
  } 
#endif", Library="DoublePendulum");
  end fmiSetInteger;
end Functions;
