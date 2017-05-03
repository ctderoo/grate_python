grate_python: a Pythonic interface to the PCGrate SX 6.1 software
__________________________________________________________________

Includes functions to initialize calculations, run scanning calculations, adapt 
model parameters in code by writing to an XML file, and extract the results from 
the output XML.

This software interface was originally written by Ryan Allured, and was adapted by 
Casey DeRoo into a python package during his Ph. D. work.

Included in the package is a minimal working example: a PCGrate project has 
been exported to XML (lamellar_grating_conical.xml) and a script using grate_python to 
rewrite that XMLf file, run the PCGrate calculation, and plot a limited set of results is
included in the /example directory. 

The best way to learn how to use this software is to contact the author(s)
directory at ctderoo7 (at) gmail.com. However, the Python code is fairly straightforward
if you're familiar with Python syntax.