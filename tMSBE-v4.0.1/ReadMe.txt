-Overview-
This ReadMe file gives background and usage instructions for the Maxwell Semiconductor Bloch Equations Model (TMSBE) 
used primarily , but not limited, to model Vertical External Cavity Surface-Emitting Lasers (VECSELs).

Full documentation for all functions is written to Doxygen (dox_config). This ReadMe does not go into each individual function
and implementation but rather highlights basic simulations, some key functions and examples, and a guide to
adding additional features into this suite, including proper testing procedures. This document should function primarily
as a quick start guide for new users.

-Prerequisites-
Assumed knowledge is basic C/C++ syntax as well as the underlying physics of this model. No derivations are given herein.
Previous work with quantum mechanics, semiconductor physics, electromagnetic waves, and geometrical optics is suggested.

-Significant Developments-
Version 1.0: (????) Full modelocking model for MSBE linear VECSEL cavity.
Version 2.0: (2018) Transverse dimension added utilizing spatial Fourier transforms
Version 3.0: (2020) TwoArm structure modules and expanded directed SBE added for non-normal incidence simulations. Tested in 1D
	3.1: (2020) Non-normal angle of incidence propagation within TwoArmCavity modules added. Tested in 1D.
	3.2: (2020) Cyclic propagation through boundary transfers added. Not yet modelocked.
	3.3: (2020) Transverse non-normal incidence modified optical path lengths
	3.5: (2020) Transverse non-normal incidence phase and V-cavity modelocking
	3.6: (2020) Fixed various bugs and issues
	3.7: (2020) Fixed various bugs and issues. Used for V-cavity manuscript
	3.8: (2020) Ring cavity implementation
	4.0: (2020) Birefringent crystal implementation

-Author Aknowledgement-
Included below are the primary authors of this coding suite along with additional authors responsbile for
signficant contributions to its development:

(????-Cur.) Isak Kilen
(2018-Cur.) Samuel McLaren
(????-Cur.) Joerg Hader
(????-Cur.) Stephan W. Koch
(????-Cur.) Jerome V. Moloney

-Future Releases-
There are obviously numerous tests that need to be done with this model. Included herein are significant future
developments that will liekly be made, in no particular order:

1.Expanded directed SBE for higher order scattering (2nd Born and onward)
2.Additional dimension 2d-TMSBE implementation
3.Variable band structures (only 2-band parabolic currently allowable)


-Basic (Non-exhaustive) Checklist for Operation-
1.Setup file
	a.Reflection computation or full modelocking?
	b.Various SBE level definitions
	c.Simulation time
	d.Output/Save frequency and duration
	e.Cavity lengths	 
2.Job/Make file
	a.Number of threads and nodes
	b.Source code location	
3.Main file
	a.Two boundaries
	b.Proper starting side
	c.Cavity elements
4.Material folder
	a.QW and Absorber material files
		aa.Density
		ab.Effectivity multiplier
		ac.Focus
	b.Cavity config
		ba.Seed pulse strength/delay
		bb.Angle of incidence

-Quickstart guide for operation-
1. 1D Single slab quantum well test
	Simplest operation is using a single slab structure with just a single refractive index and single pass over the centralized quantum well.
	Construct a single refrative index slab with a central quantum well and non-reflecting boundaries.
	a. Reflection computation holds occupation numbers constant. Run with variable seed pulse strength. What effect does this have?
	b. How does gain and GDD change with increasing pump intensities (carrier densities)
	c. Repeat with angle of incidence. Can we draw ana analogy between density and angle of incidence?
	d. Coulomb matrix is linearly interpolated according to a threshold. Show that the current state is near the converged result.
	e. What differences are there when including second Born scattering? Hint: Check temperatures as well as gain/absorption, GDD.
	f. Introduce partially reflecting boundaries to create a standing wave. The wavelength can be tuned to allow multiple modes within the cavity.
	   Form one, two, and three color cavities by altering the cavity length, central frequency, output coupling, carrier density,
2. 1D Distributed Bragg Reflector (DBR) test
	The DBR is responsible for the reflectance within the cavity. Construct a reflecting cavity with a length of air and a DBR between two absorbing
	boundaries. Initialize with seed pulses of variable strengths.
	a. Visualize the DBR stop band as a function of the relative difference between the DBR refractive indices. What about absolute differences?
	b. Is the DBR reflectance a function of the cavity length?
	c. How does DBR reflectance change with pair numbers? Show some convergence.
	d. Compare results to MATLABs multidiel function
3. Optimize AR coatings
	The DBR and QW gives a nonzero group delay dispersion (GDD) to the cavity. This is reduced using coatings on the gain chip and SESAM.
	The thickness, number of layers, and refractive indices can be optimized to flatten, center, and zero at center the GDD.
	Build a cavity in MATLAB using multidiel and optmize over 0-2 layers plus a cap layer. COmpare this result to a similar cavity in this suite.
4. Modelocking in 1D
	a. Ratio of saturable and unsaturable losses
	b. Quantum well densities
	c. Cavity length and recovery rates
5. Modelocking in a V-cavity
	a. Density tests
	b. Angle of incidence tests
	c. Relative and absolute arm lengths
6 To Be Continued...

-Test/Implementation Procedure for future release features-
0. All implementations should be tested AS THEY ARE IMPLEMENTED and should be tested in THEIR SIMPLEST FORM.
1. Standalone test comparison to analytic/experimental/confirmed reference.
2. Reflection spectrum gain, GDD comparison tests at appropriate densities
3. Modelocking simulations with convergence and instantaneous output plots for intensity, energy, power, FWHM, etc.
4. Destabilization and breaking criteria tests that highlight range of validity and required elements for proper usage.
5. Transverse operation at key parameter space points for comparisons.


Copyright (c) [2021] [Samuel A. McLaren]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
