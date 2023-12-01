# Terahertz-pumping

## Simulation File Description In Example

Completion of this example requires:
NAMD 2.10 or later (http://ks.uiuc.edu/Research/namd); VMD 1.9.3 or later (http://ks.uiuc.edu/Research/vmd; using the latest version is recommended).

1) mol.pdb: This file contains initial position coordinate data.
    
2) mol.psf: This file describes the molecular system structure.
    
par_nanochannel.prm: This file contains the parameters necessary for the molecular system to be simulated. The left-right wetting ratio of the nanochannel can be adjusted by setting the following parameters:
    
    CA     0.000000  -0.070000     1.992400
    CB     0.000000  -7.000000     1.992400
    
    
3) fix.ref: The atoms with the beta value of 1 in this file is held fixed during the simulation.
    
4) run.conf: This file is the configuration file for running computational simulations. Running the calculation. The script for simulating the cosine transform electric field has been  added to this file, with the specific code and explanation as follows:

##########################################################

tclBC		on

tclBCScript {

proc calcforces {step unique A } {

global A 

set Pi 3.1415926

###Set cosine electric field Et###

set Et [expr $A*cos(0.002*$Pi*$*$step)]

####Excludes the carbon atom from future iterations on this processor###

while {[nextatom]} { 

set chge [getcharge]

if { $chge == 0.0 } {

dropatom

continue

}

###Apply the E(t) only to the water inside the nanochannel (35>=z>=-35 angstrom)##

set var [getcoord]

set z_coord  [lindex $var 2 ];  # get the z-coordinate of the atom

if {$z_coord>=-35 && $z_coord <=35 } {

addforce [vecscale {0 1 0} [expr $Et*$chge]]};  #along the y direction

}

}

}

tclBCArgs {4.6 27000};  #A=2 V/nm, =27 THz

#############################################################

Please refer to the Namd manual for the specific meaning of the tclBCScript parameter in the following code.

5) out.dcd: This file stores the trajectory of all atom position coordinates.

## Trajectory analysis code

flow.tcl: This is a code program that calls VMD software to analyze water flow, with the command "source com_flow.tcl".

## Calculation of Water Molecule Infrared Spectrum

The IR Spectral Density Calculator Plugin  in VMD software can also be used to compute IR spectral densities.

## Reference

Qi-Lin Zhang, Tong Zhou, Chao Chang, Shi-Yu Gu, Yun-Jie Wang, Qi Liu, and Zhi Zhu, Physical Review Letters, Under review, ID: LJ18210
