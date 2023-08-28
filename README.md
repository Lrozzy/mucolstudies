# mucolstudies

Collection of scripts for performing Muon Collider studies. These scripts assume you're inside the singularity image used in the [Fermilab Muon Collider tutorial](https://mcdwiki.docs.cern.ch/tutorials/fermilab2022/computing_setup/) which does some path mapping for you on the Snowmass cluster. In particular, you should run:

`singularity run -B /collab/project/snowmass21/data/muonc:/data -B /work/$USER /cvmfs/unpacked.cern.ch/registry.hub.docker.com/infnpd/mucoll-ilc-framework\:1.6-centos8`

and then from inside the image, run:

`source /opt/ilcsoft/muonc/init_ilcsoft.sh`

(FYI, putting these two in one shell script and sourcing it will not work; it won't execute the second command until you exit the singularity. You can also update to a newer release, which you can find [here](https://confluence.infn.it/display/muoncollider/Software), but I ran into trouble running with 1.7 with my current setup.)

For makeMuonPlots.py and bad_res.py, you will need to use a different container:
Create an apptainer in your desired directory:
(NB: CARE, step 1 takes a while but only need to do it once)

`apptainer build k4toroid.sif docker://madbaron/k4test-ubuntu:latest` 

Run the apptainer in your desired directory (e.g /work/$USER instead of /home)

`apptainer run --no-home -B /collab/project/snowmass21/data/muonc:/data -B /home/$USER k4toroid.sif`

Source setup script

`source /setup.sh`


For a very simple version of reading an slcio file, see `makeMuonPlots_simple.py`. 

For a slightly more advanced version that handles making multiple histograms more elegantly, see `makeMuonPlots.py`. It may be helpful to look at these side-by-side if you're trying to learn what's going on here.

Both of these scripts will put plots in a `plots/` directory - make one if you don't have it already or it will crash.

For the SLCIO files, to get a list of the colletions you can access, run `anajob <filename>` on one of your files.
Use the `COLLECTION NAME` to access them. 
To understand what kinds of functions you can use on the particles in these collections, look at the `COLLECTION TYPE` and look up its functions [here](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/namespaceEVENT.html).

To understand in more depth how reconstruction works, you can look at the code [here](https://github.com/MuonColliderSoft/DDMarlinPandora/tree/master/src). To understand what options were passed to this code, you'll need to look at the steering files that were used to run it. Fede's are [here](https://github.com/madbaron/SteeringMacros/tree/master/Reco), but you'd have to know which one he used. For more generic ones, have a look at the ones used in the
[tutorial](https://github.com/MuonColliderSoft/MuC-Tutorial/tree/master/reconstruction).
=======
Muon Collider scripts on Snowmass
>>>>>>> fdd9f6f6efc8f3022fd87cf0cdec09c3660fc1bf
=======
Muon Collider scripts on Snowmass
>>>>>>> a08359c84050bee085be35b7fab37a7615af40a7
