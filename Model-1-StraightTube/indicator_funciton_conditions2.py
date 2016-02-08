from dolfin import *
from numpy import where

# take in wss and bdry instead of wss-bdry
# make a really simple function that takes k func and returns bdry

#def ind_func_expr(WSS_bdry):

def ind_func(bdry, WSS, interesting_domain):
	# take the WSS_bdry function and make it into vector
	# atm of getting it from the attributes in in DPspace and thats fine because so
	# is our ind function supposed to be

	wss_ = WSS.vector().array()

        # Project bdry to the wall shear stress function space to ensure
        # that both of them live in the same space. No loss of information
        # occurs here.
        bdry = interpolate(bdry, WSS.function_space())
	bdry_ = bdry.vector().array()


        assert len(wss_) == len(bdry_)

        bdry_idx = where(bdry_ > DOLFIN_EPS)[0]

        # Plotting wss_bdry is a pain, but iterpolation to DG0 seems to work
        # okish
        DG0 = FunctionSpace(WSS.function_space().mesh(), "DG", 0)
        #plot(interpolate(wss_bdry, DG0), interactive=True)

        wss_bdry = wss_[bdry_idx]
        print "Max / min WSS at boundary: {} / {}".format(max(abs(wss_bdry)),
                min(abs(wss_bdry)))

	# 50% off: thresh_L = 0.025; thresh_H = 0.075
	# 60% off:
	thresh_L = 0.000001; thresh_H = 0.04

        growth = where((abs(wss_) > thresh_H) & (bdry_ > DOLFIN_EPS))[0]
        shrink = where((abs(wss_) < thresh_L) & (bdry_ > DOLFIN_EPS))[0]

        indicator_function = Function(WSS.function_space())
        indicator_function.vector()[growth] = -1
        indicator_function.vector()[shrink] = 1

        #plot(interpolate(indicator_function, DG0), interactive=True)
        #plot(interpolate(bdry, DG0), interactive=True)

        return indicator_function
