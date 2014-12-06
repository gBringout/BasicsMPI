# BasicsMPI #

The basic functions to:
+ simulate the creation of MPI-signal using a Langevin model for the magnetization of the particle,
+ simulate the encoding princples of the MPI-signal
+ simulate a 2D FFP scanner which is made of ideal/pure fields and obtain the system matrix basis functions
+ simulate a 2D FFL scanner which is made of ideal/pure field and obtain the systemmatrix basis functions

More details on the workflow for the similation of MPI-signal for scanners is given in the [wiki](https://github.com/gBringout/BasicsMPI/wiki)

Values and ideas comes mainly from:
+ [S. Biederer, "Entwicklung eines Spektrometers zur Analyse superparamagnetischer Eisenoxid-Nanopartikel für Magnetic-Particle-Imaging", Springer, 2012](http://www.springer.com/springer+vieweg/it+%26+informatik/wissenschaften/book/978-3-8348-2406-6)
+ [T. knopp and T. M. Buzug "Magnetic Particle Imaging", Springer, 2013](http://www.springer.com/medicine/radiology/book/978-3-642-04198-3)
+ [J. Weizenecker, J. Borgert and B. Gleich, "A simulation study on the resolution and sensitivity of magnetic particle imaging",Physics in Medicine and Biology,2007](http://dx.doi.org/10.1088/0031-9155/52/21/001)
+ [J. Rahmer, J. Weizenecker, B. Gleichand and J. Borgert, "Signal encoding in magnetic particle imaging: properties of the system function",BMC Medical Imaging,2009](http://dx.doi.org/10.1186/1471-2342-9-4)
+ [G. Bringout and T. M. Buzug , "A robust and compact representation for magnetic fields in magnetic particle imaging",Biomedical Engineering / Biomedizinische Technik,2014](http://dx.doi.org/10.1515/bmt-2014-5009)

## MPI-signal creation ##
We create these graph and the associated data about the signal generation, 'in' the FFP, with the script "SignalGeneration.m"
![Alt text](/pictures/SignalGeneration.jpg?raw=true "The signal generation, 'in' the FFP")

## MPI-signal encoding ##
We create these graph and the associated data about the signal encoding, when a field offset is present, with the script "SignalEncoding.m"
![Alt text](/pictures/SignalEncoding.jpg?raw=true "The signal encoding, when an field offset is present")

## Ideal FFP scanner ##
We create these graph and the associated data about the signal generated by an 2D FFp scanner, with the script "Reco2D_IdealFFP.m"
![Alt text](/pictures/FFP_SM.jpg?raw=true "The basis functions stores in the SM of a 2D FFP scanner")

## Ideal FFL scanner ##
We create these graph and the associated data about the signal generated by an 2D FFL scanner, with the script "Reco2D_IdealFFL.m"
![Alt text](/pictures/FFL_SM.jpg?raw=true "The basis functions stores in the SM of a 2D FFL scanner")
