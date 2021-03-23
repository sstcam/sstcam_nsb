This folder contains example scripts using NSB for different investigations.

Installation instructions for joint ctapipe/nsb environment:

<<<<<<< HEAD
conda create -n nsbenv python=3.8
conda activate nsbenv
conda install pip
conda install -c cta-observatory ctapipe
pip install nsb

Configuration information:

nsb needs a mypycat.txt saved in the ~/.nsb folder to operate with the mypycat() function, use the version in this git folder. 
=======
conda create -n nsbenv python=3.8  
conda activate nsbenv  
conda install pip  
conda install -c cta-observatory ctapipe  
pip install nsb  
>>>>>>> d736a6be3b0e18cbfbd5aea7dd85dc30322abf55
