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