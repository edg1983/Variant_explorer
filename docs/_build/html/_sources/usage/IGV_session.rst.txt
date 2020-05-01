IGV session file
================
This section explains how to generate an IGV session file in the app and use it to inspect your variants
This option is active from teh "Gene details" section to inspect one or more of the candidate variants

Get the IGV session file
++++++++++++++++++++++++

From the "gene details" tab in the app select the variant(s) that you want to inspect. 
You can select any number of variants from either comound hets table and/or single variants table,
as long as they are on the same chromosome.
When you select a variant you should see the region that will be configured in IGV in a text box 
on the right of the download button.

1. Selct variant(s) of interest
2. Click on the "Download IGV session" button
3. An XML file is automatically downloaded to your computer.

Prepare for visualization
+++++++++++++++++++++++++

To visualize the variant you need to run IGV from the server, either rescomp or humbug

1. Copy the XML file downloaded from the app to any location in the server
2. Log in into the server using one of the following. Note that humbug require you are connected to VPN.

``ssh -X -l username rescomp1.well.ox.ac.uk``
``ssh -X -l username humbug.well.ox.ac.uk``

3. Load IGV using

``module load IGV/2.8.0-Java-11``

4. Launch IGV from the server

``igv.sh``

5. An IGV windows should open automatically on your computer. 
Use File > Open Session and load the XML file you have previously copied 

The IGV screen would be set automatically to the proper location, diplaying the BAM files from the selected case
labelled with affected / unaffected. Family VCF file and SV VCF will be loaded as well.
