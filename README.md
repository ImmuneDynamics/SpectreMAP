# SpectreMAP

A computational toolkit to facilitate the spatial analysis of Imaging Mass Cytometry (IMC) data in R.

![Spatial](https://wiki.centenary.org.au/download/attachments/172228252/image2021-2-25_22-32-15.png?version=1&modificationDate=1614252735692&api=v2)

Along with flow, spectral, or mass cytometry data, Spectre enables spatial analysis of Imaging Mass Cytometry (IMC) data. The Hyperion from Fluidigm (an Imaging Mass Cytometer, IMC) consists of a CyTOF (Helios) instrument, with an imaging module attached to the front. Within the imaging module, a pulsed laser scans and ablates the tissue section in incremental 1 um shots, which are then rastered together into an image, consisting of 30-40 metal signals representing different cellular or tissue markers. In order to analyse this imaging data in R, we developed an extension of Spectre, termed ‘SpectreMAP’, to import, manage, and visualise TIFF files using RStudio.

Instructions and protocols can be found at https://immunedynamics.github.io/spectre/spatial/.
