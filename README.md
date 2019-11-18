![Schmetical Diagram](https://github.com/wangshaoyun/Bond_Fluctuation_Method_SPME/blob/master/Fig1.jpg "Simulation System")
# Bond_Fluctuation_Method_SPME
## About this program
1. This is a program used to simulate weak polyelectrolytes(PE) brushes with salts by Kremer's bond fluctuation method [1,4].
2. PE brushes can be linear or star. 
3. 8*8 chains are grafted on a plate, and each chain includes 100 monomers. Same quantities counterions of PE brushes are added in the system. In addition, trivalent salts with same charge quantities as charged monomers and corresponding couterions of trivalent salts are added in the system . The total charged particles are 21333. 
3. Interaction only includes coulomb potential and Metropolis algorithm is used to update position.
4. Ewald sum is only suit for particles no more than 3000, so Monte Carlo multiple time step algorithm [2] which separt the Coulomb potential into short part and long part are used to accelerate calculate the Coulomb potential.
5. Each Monte Carlo step update the long part potential, and the long part potential which is Fourier part Coulomb potential in fact is calculated by smooth particle mesh Ewald (SPME) [3] algorithm with complexity O(Nlog(N)).
6. Multiple time step algorithm is veried with traditional Ewald Metropolis algorithm in small system with partcile 3000.
7. Constant pH titration algorithm is used to simulate the titration process.
8. Short part potential is calculated by pure cell lists algorithm with doubly linked lists.
>[1] David P. Landau, Kurt Binder. "An Guide to Monte-Carlo Simulation in Stastical Physics." 3rd ed, Cambridge Unversity Press, Cambridge UK, 2009.  
>[2] K. Bernacki, B. Hetenyi, B. J. Berne. "Multiple "Time Step" Monte Carlo Simulations: Application to charged systems with Ewald Summation." *Journal of Chemical Physics*, **121** (1), pp.44-50, 2004.  
>[3] Ulrich Essmann, Lalith Perera, Max L. Berkowitz, Tom Darden, Hsing Lee, Lee G Pedersen. "A Smooth Particle Mesh Ewald Method." *Journal of Chemical Physics*, **103** (19), pp. 8577-8593, 1995.  
>[4] I. Carmesin, Kurt. Kremer. "The Bond Fluctuation Method: A New Effective Algorithm for the Dynamics of Polym
ers in All Spatial Dimensions.", *Macromolecules*, **21** (9), pp. 2819-2823, 1988. 

## About Parameter Files 
+ energy_data.txt: It is a parameter file which includes Bjerrum length, external electric field
+ system_data.txt: It is a parameter file which includes system, time, brushes parameter.  

## Compiling files
1. Determine the parameters such as box size, brush chians, time step and so on.
2. Set the parameters in energy_data.txt and system_data.txt
3. Open terminal and enter the current file content
4. To compile the files:
```
$ make
```
  
5. To execute the files:
```
$ ./main &
```







