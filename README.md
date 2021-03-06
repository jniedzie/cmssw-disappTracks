## cmssw-disappTracks

This repo contains additions to CMSSW for the disappearing tracks analysis:

* CharginoAnalyzer - an analyzer that extracts XYZ positions of the hits and provides some additional information at different processing levels (gen, sim, reco...),

* Chargino events generation - scripts, cards, modules required to generate chargino events containing the soft pion and all information necessary for a looper reco development/analysis.

1. **Setup**

* This was tested with CMSSW_9_4_6_patch1. Go to your working directory and initialize CMSSW:

```
cmsrel CMSSW_9_4_6_patch1
cd CMSSW_9_4_6_patch1/src
cmsenv
```

* Clone this repo into CMSSW source directory _CMSSW_9_4_6_patch1/src_ (dot at the end important!):

```
git clone git@github.com:jniedzie/cmssw-disappTracks.git .
```

* Build CMSSW:

```
scram build
```

2. **Events generation: GEN-SIM step**

* Chargino events will be generated applying a MET filter of 200 GeV. Therefore, be prepared to get ~5% of the number of events specified in the command below.

* Run this scripts with proper arguments and go get a coffee (a big one if you scheduled many jobs):

```
./submitJobs.sh \
-j number_of_jobs \
-e number_of_events_per_job \
-s GEN-SIM \
-o /eos/cms/store/group/something/something_else/susy/GEN-SIM/
```

* To check if jobs are running, use `condor_q` command.

* If everything goes fine, after many hours you should see files like `chargino300GeV_ctau10cm_GEN-SIM_0.root` containing generated events in the output directory specified.

4. **Events processing: GEN-SIM-RAW, GEN-SIM-RAW-RECO and miniAOD steps**

* Once GEN-SIM events are ready, you can run next step using the same script, but this time specifying input path. No need to specify number of events anymore - all events in the file will be processed:

```
./submitJobs.sh \
-j number_of_jobs \ 
-s GEN-SIM-RAW \
-i /eos/cms/store/group/something/something_else/susy/GEN-SIM/ \
-o /eos/cms/store/group/something/something_else/susy/GEN-SIM-RAW/
```

In order to perform further steps, just replace GEN-SIM-RAW with GEN-SIM-RAW-RECO or miniAOD (and remember to update input and output paths).

5. **Some random info that may be useful**

* make sure that you are in a BASH shell,
* if you have some problems with jobs submission, try adding something like this to your ~/.bashrc file:
`export X509_USER_PROXY=/afs/cern.ch/user/a/aalibaba/x509up_u12345`
changing the path to point to an existing file in your user directory.

