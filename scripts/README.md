This folder contains sandbox scripts using nsb for different investigations.

Very hacky installation instructions for joint ctapipe/nsb/photutils environment:

conda create -n nsbenv python=3.8  
conda activate nsbenv  
conda install pip  
conda install -c cta-observatory ctapipe  
pip install nsb  
conda install photutils  
pip uninstall numpy  
pip install numpy  
pip install starfield (for availability plots).

Configuration information:

nsb needs a mypycat.txt saved in the ~/.nsb folder to operate with the mypycat() function, use the version in this git folder. 

Usage Instructions:

Generate a gaia catalogue using the instructions on the nsb pip webpage, note this must be done interactively.

Add RA/DEC of desired source to ~/.nsb/mypycat.txt.

Use skycoord.py to obtain RA/DEC/ALT/AZ in correct format for other scripts. Time must be provided in UTC.

Use savefield.py to write observation fits files (with a runname argument parser argument, i.e. python savefield.py test1) using nsb. Note parameters google doc file here: https://docs.google.com/spreadsheets/d/14i0Uk78AZhnojA8j56SqrpoPtfSDufpp_r0ehV3PyO0/edit?usp=sharing. Set observation time, alt, az and source at the top of the script, and an nsb config file will be automatically generated.

Use fromfits.py to perform aperture photometry on the generate fov maps from savefield.py, assuming a CHEC-like camera architecture (some parameters need updating), this generates fov and pixel plots, as well as an astropy table file to be used with the other scripts. Ensure observation time and RA/DEC are present and correct.

traybake.py and difference.py are basic scripts to generate plots per superpixel and TM, difference.py takes two astropy tables and computes the difference per superpixel and TM.

timespan.py is a script to generate nsb timescale plots, in both nLb and Hz given certain assumptions.

obstime.py is a script to calculate the possible observing time gain for a particular source, given an nsb config file (though that doesn't do very much given how nsb is set up).