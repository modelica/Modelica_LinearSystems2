within Modelica_LinearSystems2.WorkInProgress.Tests.Internal;
function clock

input Real dummy=0;
output Real runtime_ms;

external "C" runtime_ms = clock_ms(dummy);

annotation(Include = "
#include <time.h>
double clock_ms(double dummy)
{
  return (double)clock() * 1000 / (double)CLOCKS_PER_SEC; 
 
}
 
");

end clock;
