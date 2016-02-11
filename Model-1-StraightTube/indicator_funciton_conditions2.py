from dolfin import *
from numpy import where

# take in wss and bdry instead of wss-bdry
# make a really simple function that takes k func and returns bdry

#def ind_func_expr(WSS_bdry):

def ind_func(bdry, WSS, interesting_domain):
	
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

	thresh_L = 0.000001; thresh_H = 0.04

        growth = where((abs(wss_) > thresh_H) & (bdry_ > DOLFIN_EPS))[0]
        shrink = where((abs(wss_) < thresh_L) & (bdry_ > DOLFIN_EPS))[0]

        indicator_function = Function(WSS.function_space())
        indicator_function.vector()[growth] = 1
        indicator_function.vector()[shrink] = -1

	#from IPython import embed; embed()

	# Naar interpolerer sa tar man verdien som er midt pa cellen fx 1/3
	# Magne tegnet det 
	# careful about interpolating from DG0 to DG1 
	

        #plot(interpolate(indicator_function, DG0), interactive=True, title = 'indicator function in DG0')
        #plot(interpolate(bdry, DG0), interactive=True)
	


        return indicator_function
