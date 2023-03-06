# AMBER-ALPIDE 2022 Correlator
Framework used to correlate AMBER DAQ and ALPIDE DAQ
during Q4 2022 CERN test experiment.

**Alpide side** Supply a raw timestitched file.
**AMBER side** Take a trlo ROOT file that contains usual `runNumber`, `spillNumber`, `eventNumber`, `eventTime` branches.
The AMBER ROOT file entries have to be sorted first with ascending eventNumber (important!). To do this just pass it to the ``sortByEventNumber.C`` script by editing the `string fName` variable (first line in macro function).

Usage:
```sh
./main --amber=FILE_AMBER.root --alpide=FILE_ALPIDE.root
```

C++17 could be required to compile the target.
To check C++ version currently underlying your ROOT execute the following:
``which root-config | xargs cat | grep cxxversionflag=``
Back-up the root-config shell script and change the version to 17.
