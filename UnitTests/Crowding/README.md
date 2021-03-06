# Crowding influenced bimolecular association

In a crowded environment often macromolecular species influence other proteins' reactions indirectly, e.g., through volume exclusion and diffusion inhibition. On the one side, the crowding particles occupy volume and, thereby, reduce the volume available to the reactants which leads to an increase of the reaction rates. On the other side, crowding slows down the diffusion of particles and, thereby, decreases the rate of collision, and consequently the reaction rates. The crowding benchmark has been suggested in at a Dagstuhl seminar (Andrews et al. 2015) to test various spatial simulators. Already in these first tests it has been shown that the simulation results vary significantly. All simulations showed for the (nearly) diffusion limited case that the reaction rate decreases with increasing crowding. In the activation limited case, no unique qualitative behavior pattern can be derived. The results vary reflecting the idiosyncrasies of the respective collision handling within the implemented methods. In our calculations the monotonic decrease in the diffusion-limited case could be reproduced. The quantitative results again show some deviation between the different methods which depends on how the excluded volumes are realized.  Due to the lack of an analytic solution in this cases and thus of a ground truth, we expect all methods to show with increasing crowding a qualitative decrease of the reaction rate.    

All simulations are executed in a box of 50 * 50 * 50 nm with periodic boundary conditions. The solely reaction is A+B->B+C, all particles have a radius of 0.5 nm, and a diffusion coefficient of 10 nm²/µs. The simulation starts with 1000 A and B, and a number of C that is calculated depending on the crowding degree. Please note that depending on the realization of excluded volumes the max. degree of crowding that can be realized varies. For example in Smoldyn experiments have been run with close to 100% crowding (see Dagstuhl report), whereas based on non-overlapping spheres the highest density is around 74%. In our approach, we only considered crowding of up to 20% in 5% steps. 10 replications were run. In the nearly diffusion limited case each collision resulted in a reaction, whereas for the activation limited case we assumed that 20% of collisions lead to a reaction.  

Experiment setting for eGFRD (scaled down): 

All simulations were performed in a box of the size (23.21nm)^3 = (50nm * 10^(-1/3))^3.
The particles A, B and C have a radius of 0.5 nm and a diffusion coefficient of 10^-11 m^2/s. The first reaction rule is

A+B -> A-B @ ki = 85*10^-21 m^3/s

where ki is the intrinsic rate. A and B form the complex A-B (r = 0.5 nm, D = 0), because eGFRD only allows one product in a bimolecular reaction. Hence we need a second reaction rule.

A-B -> B + C @ k = inf

Since the rate k is infinity (those are possible for unary reactions), the complex A-B decays after a short time (it looks like A-B exists only for a few ns at most). The simulation starts with 100 A and 100 B particles. The number of C particles is chosen in a way, that A+B+C fill x % of the box (crowding factor). 
