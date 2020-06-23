Installation
============

Download the app
++++++++++++++++

Download the latest release of the app from github

Access `release page <https://github.com/edg1983/Variant_explorer/releases>`_ and download the latest zip archive

Extract the archive in any folder you want. 
A new folder named Variant_explorer-XXX will be created automatically

Prepare your system
+++++++++++++++++++

Windows system
--------------

Download and install the following programs if you don’t have them already

1. `R v3.6.3 <https://cran.rstudio.com/bin/windows/base/R-3.6.3-win.exe>`_
2. `Rtools v3.5 <https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe>`_
    When prompted for installation path use a folder name without blank spaces
    Check the option “Add Rtools to PATH” when asked
3. `R studio v1.2 win <https://download1.rstudio.org/desktop/windows/RStudio-1.2.5033.exe>`_

MacOS system
------------

Download and install the following programs if you don’t have them already

1. R v3.6.3
    `Normal version <https://cran.rstudio.com/bin/macosx/R-3.6.3.nn.pkg>`_ or
    `Catalina version <https://cran.rstudio.com/bin/macosx/R-3.6.3.pkg>`_
2. `R studio v1.2 mac <https://download1.rstudio.org/desktop/macos/RStudio-1.2.5033.dmg>`_

Linux system
------------

Ensure you have the following lib installed

- libcurl4-openssl-dev 
- libssl-dev 
- libsodium-dev
	
You can install using

``sudo apt-get install libcurl4-openssl-dev libssl-dev libsodium-dev``

Run the app for the first time
++++++++++++++++++++++++++++++

1. When you have installed the required software, open R studio. 
    In the files panel at the bottom right navigate to the folder where you have installed the app. 
    Click more > Set as working directory
2. Open the app.R file by clicking on it
    Now in the top-left panel, you will see the app source code.
3. Click on the “Run app” button.
    Answer yes if prompted to install new shiny packages or compile from source.

The first time you run the app, it will install all the additional packages needed, so it could require some time.
When the installation is done, the app will start automatically