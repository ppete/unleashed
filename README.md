# The Unleashed Concurrency Library (UCL)

* ./include:     contains the UCL library files.
* ./scripts:     contains shell scripts that compile test UCL code on various systems
* ./examples:    contains test and benchmark codes for the UCL framework
    * ./android:   sample AndroidStudio project that uses UCL. Note AndroidStudio 1.3 and gradle-2.5 are required for testing the simple application.
    * ./tasks:     test applications using the UCL-task framework
    * ./container: tests for several containers (skiplist, list, queue, stack) and memory management approaches (for both HTM and non-HTM systems), including garbage collection, publish and scan techniques, epochs, ..
    * ./locks:     sample codes using UCL's lock and transactional elision implementations
* ./util:        utility codes


### Integration on Android - UCLTest.tar.bz2

For installation, AndroidStudio 1.3, gradle-2.5 are required. 
When importing the project, the location of gradle-2.5 and the SDK may need to 
be modified. The project expects a environment variable UCL_HOME be set to 
UCL's top directory. For running the code on an actual device, the device 
needs to support NDK >= r10, otherwise some standard library functions cannot be 
found (e.g., rand).

### Contributors
* Nicholas Dzugan
* Christina Peterson
* Peter Pirkelbauer
* Amalee Wilson
