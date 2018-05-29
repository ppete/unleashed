The BLAZE Concurrent Data Structure Library
-------------------------------------------

./include: contains the BLAZE library files.

./example: contains sample code and a sample AndroidStudio project that uses
           BLAZE (BlazeTest.tar.bz2). 
           
./scripts: contains test scripts.

Integration on Android - BlazeTest.tar.bz2
------------------------------------------
For installation, AndroidStudio 1.3, gradle-2.5 are required. 
When importing the project, the location of gradle-2.5 and the SDK may need to 
be modified. The project expects a environment variable BLAZE_HOME be set to 
BLAZE's top directory. For running the code on an actual device, the device 
needs to support NDK >= r10, otherwise some standard library functions cannot be 
found (e.g., rand).

Contributors:

Nicholas Dzugan
Christina Peterson
Amalee Wilson
