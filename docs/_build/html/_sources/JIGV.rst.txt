Use the JIGV script
===================
This section explain how to use JIGV to visualize selected variants

This procedure allows you to load an IGV-like view of the data on your browser, without having to copy BAM/VCF files.
This function use jigv program to create a visualization service on the server. 
Using connection tunneling you can then link this service to your local machine to see data directly in your browser

Configure Putty
+++++++++++++++
This need to be done just once

Follow the guide in `this webpage <https://www.skyverge.com/blog/how-to-set-up-an-ssh-tunnel-with-putty/>`_

Configuration parameters for putty:

- Source port: 5001
- Destination: localhost:5001

Generate the script
+++++++++++++++++++

From the "gene details" tab in the app select the variant(s) that you want to inspect and then click the "Download JIGV scrip" button

A bash (.sh) script is generated automatically and downloaded to your computer.

Prepare for visualization
+++++++++++++++++++++++++

On server side
--------------

1. Access the server (eg. humbug) and copy the downloaded script to any location
2. Run the script using: ``bash downloaded_script.sh``
3. Leave the server session open and proceed on your computer

On your computer
----------------
Windows
~~~~~~~
1. Open Putty and activate the previously configured connection
2. Leave putty active
3. On your browser navigate to localhost:5001

Linux / Mac
~~~~~~~~~~~
1. Open the terminal and use: ``ssh -N -L localhost:5001:localhost:5001 username@humbug.well.ox.ac.uk``
2. Leave the terminal open
3. On your browser navigate to localhost:5001
