# APEX_G4MC
This is a modified version of G4MC, where the magnetic fields are given according to SNAKE package, fringe fields are also added, apertures are updated (from SIMC), several new volumes are added to simulate a more realistic acceptances. For APEX mode simulations, septum magnet and it's vaccum chambers are also defined.


Magnetic Fields:

Only dipole central field is using the fortran code (dipolert2.f) which is a part of the SNAKE package.
For the other fields (quads central and fringe fields, dipole fringe field) SNAKE functions are translated to cpp-language and implemented in the G4-codes (a few years ago I found this way more effective, but didn't manage to translate the last part - dipole central field).
In the current version the fields are given as global field (by HRSEMField.cc).
In old versions the fields were localized in volumes (see HRSEMFieldSetup.cc). This way was not very productive when I added fringe fields and tryed to run simulations with different step sizes, and also monitor the field values on each of the track steps.


Generator:

The initial particle generation is done by "HRSPrimaryGeneratorAction.cc".
It uses geant macro input files, while for APEX full simulation, I used root files as input to generate e- e+ pairs... 
We can modify the particles generation with different ways - depends on what process (particles) you want to generate


Detector:

There are two modes to construct the detector.
1. HRS standard mode - where you have vacuum, Q1-quad, Q2 quad, Dipole, Q3-quad
2. APEX mode - where we have displaced targets (by about dz=-100 cm), septum magnet with it's vacuum chambers and then the standard HRS spectrometers (Q1, Q2, Dipole, Q3).
Switching the septum "ON" in Detector.ini, the septum mode is built automatically.
With "SeptumOn=0" the standard HRS can be used.
3. The Sieve slits is defined for APEX modes, but it needs one line modification to place it in the HRS standard mode place.
the following parameters (from Detector.ini) are used the most 

SetupLSieveSlit=1; <br>
SetupRSieveSlit=0; <br>
sieve_thickness=0.5; <br>
Rad_tail_foil_postn=-104.8; <br>
#if fringe fileds are not needed FringeField=0 <br>
FringeField=1; <br>
SeptumOn=1; <br>
SeptumNew=1; <br>
SeptumFieldScale=1.0; <br>
Q1Sos=1; <br>
LFocalPlaneAngle=-45.0 <br>
LHRSAngle=12.5; <br>
RHRSAngle=347.5; <br>

The central momentums are defined in HRSUsage.ini  <br>
BeamEnergy  LHRSMomentum  RHRSMomentum  BeamTiltedAngle #name list <br>
2200          1500         1500         0.0             #value list <br>
These are the only numbers I changed in HRSUsage.ini <br>
The rest I inherited from old versions and haven't changed (probably most of them are not useful anymore).


Simulation Analyses:

I don't use sensitive detectors in my analyses, instead of that I directly analyze track G4Steps (see HRSSteppingAction.cc).
I use ASCII output files to write the required information during the simulations.
The outputs easily can be changed to root trees, 
but ASCII files are better when you run several hundreds of parallel jobs and they all write information in the same file simultaneously.

There are a bunch of planes (defined in DetectorConstruction and used in SteppingAction to follow the particles) that I defined and record the track coordinates on these planes.
Like "Q1Front", "vdc1Plane", "RSvSlBack" and "LSvSlBack" right and left Sieve back planes
I haven't deleted them - may be useful.


How to run:

1. Set up geant: <br>
A locally installed geant4 version should be ok. <br>
On ifarm (@JLab) <br>
> setenv JLAB_ROOT /site/12gev_phys <br>
> source $JLAB_ROOT/softenv.csh 2.0 <br>


2. Download the package to your local directory 
> git clone https://github.com/kvardan/APEX_G4MC <br>


3. Build:
> cd APEX_G4MC <br>
> mkdir build <br>
> cd build <br>
> cmake ..  <br>
> make <br>


4.Run:
After above steps an exacutable "G4MC" will appear in build/ directory
> cp G4MC ../test/no_sieve_central_target/files_in/ <br>
> cd ../test/no_sieve_central_target/run <br>
> chmod 744 run_0003.sh <br>
> ./run_0003.sh <br>

Sometimes it it useful to see the detector visually <br>
You can go to files_in/ directory and type: 
> ./G4MC -i -x 1
