### How to Use Geant code
This code is for E73 Geant4 simulation<br>

start the run by the following commands. You can get 10000 events in this case.
```sh
$ ./bin/Linux-g++/RCSim
$ /control/execute/vis.mac
$ /run/beamOn 10000
```
### Notes
- You can change initial momentum values at src/PrimaryGeneratorAction.cc
- Run logs will be saved as rootfile. You can edit this in src/RunAction.cc


### How to use analysis code
- Run ana/method1.cc to calculate initial momentum in analytical method <br>
 ![スクリーンショット 2024-06-08 0 11 40](https://github.com/Kohki-Amemiya/E73_ReconstMom/assets/144120249/8ab206ea-30cb-4574-8f54-84f81679795a)
  
- Run ana/method2.cc to calculate initial momentum by beta scanning method <br>
![スクリーンショット 2024-06-08 0 12 30](https://github.com/Kohki-Amemiya/E73_ReconstMom/assets/144120249/b7db2630-cf70-4e35-a356-5d38edd569bf)

  
