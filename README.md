# Installation #

Source the openfoam environment (bashrc or cshrc file) to load the proper environment variables. Place this folder in the $WM_PROJECT_USER_DIR/utilities/ directory and run wmake libso from the source directory. The application will be installed in the directory set up by the environment variable $FOAM_USER_LIBBIN

# Usage #

In order to use this boundary condition make sure to include liblaminarBLInlet.so in the libs section of the controlDict file. The flowspeed(velocity) and delta(boundary layer height) must be defined. The typical usage is as follows:-

```
myPatch
{
    type         laminarBLInlet;
    flowSpeed    (3.1 0 0); // Free-Stream velocity
    delta        0.1;       // boundary layer height at the patch location
    value        uniform (0 0 0); // default initial value, must be provided
}
```