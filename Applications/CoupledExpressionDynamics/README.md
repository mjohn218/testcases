# Coupled Expression Dynamics
Oscillations in protein expression levels can exhibit sensitivity to copy number fluctuations that are only captured in stochastic solvers such as Gillespie or single-particle RD. Based on a simple model of circadian oscillations (Vilar et al, PNAS 2002), we have simulated the behavior of an activator protein A and repressor protein R that are produced from mRNA transcribed from a single copy of a gene (one for each protein). Coupling of A and R expression is driven by positive feedback of the activator A, which binds to each gene’s promoters to enhance transcription. Protein R also binds to A to effectively degrade it, and all proteins and mRNA are degraded. If the degradation rate of protein R is slow, the oscillations will end in the deterministic model, but persist in the stochastic method. 


Species Initial Concentration Diffusion Coefficient

1. A 0 uM (0 molecules) 10 um2/s
2. R 0 uM (0 molecules) 10 um2/s
3. C 0 uM (0 molecules) 10 um2/s
4. PrmA 0.000397 uM (1 molecule) 10 um2/s
5. PrmA_bound 0 uM (0 molecules) 10 um2/s
6. PrmR 0.000397 uM (1 molecule) 10 um2/s
7. PrmR_bound 0 uM (0 molecules) 10 um2/s
8. mRNA_A 0 uM (0 molecules) 10 um2/s
9. mRNA_R 0 uM (0 molecules) 10 um2/s

Reactions

1. A + R → C (1204 uM-1s-1) 
2. PrmA → PrmA + mRNA_A (50 s-1)
3. PrmA_bound → PrmA_bound + mRNA_A (500 s-1)
4. PrmR → PrmR + mRNA_R (0.01 s-1)
5. PrmR_bound → PrmR_bound + mRNA_R (50 s-1)
6. PrmA + A ↔ PrmA_bound (602 uM-1s-1) (50 s-1)
7. PrmR + A↔ PrmR_bound (602 uM-1s-1) (100 s-1)
8. mRNA_A → mRNA_A + A (50 s-1)
9. mRNA_R → mRNA_R + R (5 s-1)
10. mRNA_A → NULL (10 s-1)
11. mRNA_R → NULL (0.5 s-1)
12. A → NULL (1 s-1)
13. R→ NULL (0.2 s-1)
14. C→ R (1 s-1)

Parameters

V=4.189um^3 (Sphere with R=1um). 
Run time=200 s. 

Alternative systems:
Oscillations are eliminated in deterministic simulations when protein R degradation is slowed from 0.2 to 0.05 s-1: Rxn #13, R→ NULL (0.05 s-1)


# RESULTS:
For a Volume of 4.189um^3 (Sphere with R=1um), the ODE and PDE give identical results. The stochastic simulations (Gillespie) produce the same average trends in oscillations, results from 10 trajectories of 200s each.  
Single particle: 
NERDSS is the same for D=10. In the newer version, unimolecular reactions are treated as a population, which is more accurate than one-at-a-time for fast reactions. dt=10us.  Smoldyn is not significantly differen.
However, MCELL has distinctly shorter oscillations (of significance), but with a similar lag-time.  

# Statistics: 
Calculate the average time separation between maxima in A and R, and the lagtime between peaks in A and R.
These were calculated in two ways. 
1)Identifying peaks and then measuring the distance between them. Same for Lag-time. MATLAB PROGRAM: calc_peak_seps.m
2) Using a discrete FFT, and using the maximum coefficient to idenify the wavelength of the peaks. The cross-correlation between the A and R time-series was used to identify the lag-time. For stochastic simulations, the mode value from multiple trajectories was taken. MATLAB PROGRAM: findpeaksfftSTOCH.m --works for any simulation type.
# Table
|  | A wavelength | R wavelength | A-R lagtime (Cross-correlation) |
|---|---|---|---|
|1. ODE (1 traj, 8 peaks) | 25.2s +/- 0.1 SEM | 25.1s +/- 0.1 SEM | 6s +/-0.1 SEM |
|2. ODE & PDE FFT Mode <br> (1 traj, 8peaks, 200s+zeropad to 5000s) | 25.12s | 25.12s | 6.55s |
|1. PDE (1 traj, 8 peaks) | 25.3s +/- 0.2 SEM | 25s +/- 0 SEM | 5.9s +/-0.2 SEM |
|1. SSA (10 traj: 74 and 72 peaks) | 25s +/- 0.4 SEM | 25s +/-0.4SEM | 5.98s +/- 0.1 SEM | 
|2. SSA  (10 traj, 200s +zeropad to 5000s) FFT Mode | 24.8s | 25.1s | 6.63s | 
|1. Smoldyn  (1traj 8 peaks) | 25.9s   +/- 1.3s SEM | 25.6s +/- 1.3s SEM | 5.62s +/- 0.3 SEM | 
|2. Smoldyn (1traj, 200s + zeropad to 5000s) FFT | 26.6s | 26.3s | 6.3s | 
|1. Mcell (1traj, 9peaks) | 21.8s+/- 0.7 SEM | 21.8s +/- 0.7 SEM | 6.3s+/- 0.4 SEM| 
|2. Mcell  (1traj, 200s + zeropad to 5000s) FFT | 21.8s | 21.8s | 6.64ss | 
|1. FPR1 (1traj, 1000s+zeropad to 5000s) FFT | 24.3s | 24.3s | 6.45s | 
|2. FPR1 (1traj, 1000s 42 peaks) | 23.6s+/-0.6 (SEM) | 23.6s+/0.6s (SEM) | 5.9s+/-0.1s (SEM) |
|3. NERDSS (1traj, 600s+zeropad to 5000s) FFT | 25.2s | 25.2s | 6.6s | 
|4. NERDSS (1traj, 600s 26 peaks) | 25.4s+/-0.5 (SEM) | 25.4s+/0.6s (SEM) | 6s+/-0.1s (SEM) |


# LOCALIZED PRIMERS:
To introduce a spatial gradient in concentrations, 
# geometry:
the genes for A and R were both centered in a sphere of R=0.1um, still with a single molecule of each
# initial conc.
prmA = 1 molecule, New localized Conc.=0.3964uM
prmR = 1 molecule, New Localized Conc.=0.3964uM


# diffusion
The genes were uniform in this subvolume, and no barrier to them existed. 
prmA: D=0.0um^2/s
prmA_bound:D=0.0um^2/s
prmR: D=0.0um^2/s
prmR_bound: D=0.0um^2/s

They could not diffuse, however. 
# Sim 1:
DA=5um2/s, DR=50, DRNA/R=5, DC=4.5, D_promoters=0
Periphery
|P. PDE (1 traj, 9 peaks) | 24.5s +/- 0.3 SEM | 24.3s +/- 0.1 SEM | 5.9s +/-0.2 SEM |
Center
|C. PDE (1 traj, 9 peaks) | 24.5s +/- 0.3 SEM | 24.3s +/- 0.1 SEM | 5.9s +/-0.2 SEM |

# Sim 2:
DA=2um2/s, DR=50, DRNA/R=5, DC=1.9, D_promoters=0
Periphery
|P. PDE (1 traj, 10 peaks) | 22.6s +/- 0.3 SEM | 22.5s +/- 0.1 SEM | 5.8s +/-0.2 SEM |
Center
|C. PDE (1 traj, 10 peaks) | 22.6s +/- 0.3 SEM | 22.5s +/- 0.1 SEM | 5.8s +/-0.2 SEM |

# DECAY R set to 0.05:
In this case, no oscillations occur in a deterministic model. They occur in a stochastic sim. The significantly longer oscillation time is due to the much slower decay of the R, the spike in A is not so different. 
NERDSS:
|NERDSS FPR (1traj, 1000s 16 peaks) | 63s+/-2.8 (SEM) | 63s+/-2.8s (SEM) | 7.6+/- 0.3s (SEM) |


# results
The results from the PDE are essentialy unchanged. There is a very small decrease in the height of the A and R oscillations, but it could probably be attributed to the numerical issues associatied with the coarse mesh (~0.1um). The difference between the copies of A or R proteins in the center, versus on the periphery, are barely distinguishable, see Figure: AR_oscillations_vsTime_PDE_centeredGenes.eps

However, with slower D_A, The A copies did start to concentrate more in the center, a small but visible amount. This result was not sensitive to D_R, only D_A. R showed almsot indistinguishable differences in copy numbers throughout the cell, remaining well mixed for all models. 

Smoldyn was set up for this geometry, but it did not work!

# conclusion
The overall size of the cell (R=1um) is small enough that having to reach a specific point in the center of the cell to transcribe a gene does not affect the oscillations, unless D_A is slowed markedly. Oscillations then become faster, due to increased localization of A in the center. 

References: Vilar, J. M.G. et al, Mechanisms of noise-resistance in genetic oscillators. PNAS, 99:5988-5992 (2002).
