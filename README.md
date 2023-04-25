# AMBER-ALPIDE 2022 Correlator
Framework used to correlate AMBER DAQ and ALPIDE DAQ
during Q4 2022 CERN test experiment.

**Alpide side** 
Due to problems with our timesorter, some readout cycles were split over multiple entries. 
Take a raw stitched file (before any clustering is done) and pass it to the ``sortAlpideFile.C`` script by editing the `string fName` variable (first line in the macro function).

**AMBER side** 
Take a trlo ROOT file that contains usual `runNumber`, `spillNumber`, `eventNumber`, `eventTime` branches.
The AMBER ROOT-file entries have to be sorted first with ascending eventNumber. To do this just pass it to the ``sortByEventNumber.C`` script by editing the `string fName` variable (first line in the macro function).

Usage:
```sh
./main --amber=SORTED_FILE_AMBER.root --alpide=SORTED_FILE_ALPIDE.root
```

C++17 could be required to compile the target. If you're having compilation problems,
check C++ version currently underlying your ROOT by executing the following:

``which root-config | xargs cat | grep cxxversionflag=``

Back-up the root-config shell script and change the version to 17.
