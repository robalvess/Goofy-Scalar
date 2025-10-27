import sys
sys.path.append('/home/alves/pyrate/results/Goofy_Neutrino/PythonOutput')

from Goofy_Neutrino import RGEsolver

##############################################
# First, create an instance of the RGEsolver #
##############################################

rge = RGEsolver('rge', tmin=0, tmax=20, initialScale=0)


##########################################################
# We fix the running scheme and initial conditions below #
##########################################################

# Running scheme :

rge.loops = {'GaugeCouplings': 3,
             'Yukawas': 2,
             'QuarticTerms': 2,
             'ScalarMasses': 2,
             'Vevs': 2}

# Gauge Couplings

rge.g_U1Y.initialValue = 0
rge.g_SU2L.initialValue = 0

# Yukawa Couplings

rge.Ynu.initialValue = [[0., 0., 0.],
                        [0., 0., 0.],
                        [0., 0., 0.]]

rge.lambdaN.initialValue = [[0., 0., 0.],
                            [0., 0., 0.],
                            [0., 0., 0.]]

rge.lambdaR.initialValue = [[0., 0., 0.],
                            [0., 0., 0.],
                            [0., 0., 0.]]

rge.lambdaL.initialValue = [[0., 0., 0.],
                            [0., 0., 0.],
                            [0., 0., 0.]]


# Quartic Couplings

rge.lambdaH.initialValue = 0
rge.lambdaP.initialValue = 0
rge.lambdaS.initialValue = 0
rge.lambda1.initialValue = 0

# Scalar Mass Couplings

rge.m1.initialValue = 0

# Vacuum-expectation Values

rge.vH.initialValue = 0
rge.vS.initialValue = 0


############################
# Solve the system of RGEs #
############################

rge.solve(step = .05)

# Another way to call rge.solve() :
# rge.solve(Npoints = 500)

####################
# Plot the results #
####################

rge.plot(subPlots=True, printLoopLevel=True)


#############################################
# Possibly save the results for a later use #
#############################################

# Save the results in some file

# rge.save('rgeResults.save')

# Later, load the rge object with :

# rge = RGEsolver.load('rgeResults.save')

