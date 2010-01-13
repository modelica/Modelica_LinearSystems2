within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function clock2

output Real runtime_ms;

external "C" clock_ms(runtime_ms);

annotation(Include = "
#include <time.h>
int clock_ms(double *rt)
{
*rt=(double)clock() * 1000 / (double)CLOCKS_PER_SEC; 
 return 0;
}
 
");

end clock2;
