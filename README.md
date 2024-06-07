### How to Use Geant code

This code is for E73 Geant4 simulation

1. enter the command below first.

```sh
$ ./bin/Linux-g++/RCSim
```

2. you can also set up the GUI.

```sh
$ /control/execute/vis.mac
```

3. start the run, and you can get 10000 events.

```sh
$ /run/beamOn 10000
```

### How to use analysis code
- Run ana/method1.cc to calculate initial momentum in analytical method
  
- Run ana/method2.cc to calculate initial momentum by beta scanning method
