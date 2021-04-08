This folder contains example scripts using NSB for different investigations.

Very hacky installation instructions for joint ctapipe/nsb/photutils environment:

conda create -n nsbenv python=3.8  
conda activate nsbenv  
conda install pip  
conda install -c cta-observatory ctapipe  
pip install nsb  
conda install photutils  
pip uninstall numpy  
pip install numpy  

Configuration information:

nsb needs a mypycat.txt saved in the ~/.nsb folder to operate with the mypycat() function, use the version in this git folder. 

Usage Instructions:

Use savefield.py to write fits files (with a runname argument parser argument, i.e. python savefield.py test1) using nsb. Note parameters google doc file here: https://docs.google.com/spreadsheets/d/14i0Uk78AZhnojA8j56SqrpoPtfSDufpp_r0ehV3PyO0/edit?usp=sharing. Set observation time and source in the script, making sure the nsb file my_config.cfg matches.

Use fromfits.py to perform aperture photometry on the fov maps, assuming a CHEC-like camera architecture (some parameters need updating).