def build_dipping_plate(xmin,xmax,ymin,ymax,cellsize,xhinge,horizon_dip,z_base):
    # some nifty description about this function

    import numpy as np 
    
    # create x-y mesh
    xvec = np.linspace(xmin,xmax,int((xmax+cellsize-xmin)/cellsize))
    yvec = np.linspace(ymin,ymax,int((ymax+cellsize-ymin)/cellsize))
    XX,YY = np.meshgrid(xvec,yvec)

    # build z-surface
    ZZ = np.ones(XX.shape) * z_base
    ZZ2 = z_base - (XX-xhinge)*np.tan(np.deg2rad(horizon_dip))
    ZZ[XX>=xhinge] = ZZ2[XX>=xhinge]

    # output the file
    outputfilename = 'basement_surface_%.1f_degrees.csv' % horizon_dip
    hdr = 'x-lims %d %d, y-lims %d %d, cellsize %d, xhinge %d, dip %d, z_base %d' % (xmin,xmax,ymin,ymax,cellsize,xhinge,horizon_dip,z_base)
    np.savetxt(outputfilename,ZZ,delimiter=",",header=hdr,fmt='%.3f')

def build_effective_susc_uniform(xmin,xmax,ymin,ymax,cellsize,crust_susc,rem_inc,rem_dec,rem_susc,field_inc,field_dec):

    import numpy as np
    from discretize.utils import mkvc, refine_tree_xyz
    from SimPEG.utils import plot2Ddata, model_builder, surface2ind_topo, mat_utils

     # create x-y mesh
    xvec = np.linspace(xmin,xmax,int((xmax+cellsize-xmin)/cellsize))
    yvec = np.linspace(ymin,ymax,int((ymax+cellsize-ymin)/cellsize))
    XX,YY = np.meshgrid(xvec,yvec)

    # Compute the unit direction of the inducing field in Cartesian coordinates
    field_dir = mat_utils.dip_azimuth2cartesian(field_inc, field_dec)

    # Multiply susceptibility model to obtain the x, y, z components of the
    # effective susceptibility contribution from induced magnetization.
    susc_model = np.outer(crust_susc, field_dir)

    rem_model = (rem_susc * mat_utils.dip_azimuth2cartesian(rem_inc, rem_dec))
    
    # Define effective susceptibility model as a vector np.r_[chi_x, chi_y, chi_z]
    plotting_model = susc_model + rem_model
    plotting_model = np.sqrt(np.sum(plotting_model, axis=1) ** 2)

    # populate zmesh
    conversion_factor = 79617.8 # convert to ucgs
    ZZ = np.ones(XX.shape) * plotting_model * conversion_factor

    # output the file
    outputfilename = 'basement_susceptibility_%.1f.csv' % rem_susc
    hdr = 'x-lims %d %d, y-lims %d %d, cellsize %d, crust-susceptibility %.2f, field-inc %.2f, field-dec %.2f, rem-susceptibility %.2f, rem-inc %.2f, rem-dec %.2f' % (xmin,xmax,ymin,ymax,cellsize,crust_susc,field_inc,field_dec,rem_susc,rem_inc,rem_dec)
    np.savetxt(outputfilename,ZZ,delimiter=",",header=hdr,fmt='%.3f')

def build_effective_susc_stripe(xmin,xmax,ymin,ymax,stripe_width,cellsize,crust_susc,rem_inc,rem_dec,rem_susc,field_inc,field_dec):

    import numpy as np
    from discretize.utils import mkvc, refine_tree_xyz
    from SimPEG.utils import plot2Ddata, model_builder, surface2ind_topo, mat_utils

     # create x-y mesh
    xvec = np.linspace(xmin,xmax,int((xmax+cellsize-xmin)/cellsize))
    yvec = np.linspace(ymin,ymax,int((ymax+cellsize-ymin)/cellsize))
    XX,YY = np.meshgrid(xvec,yvec)

    # Compute the unit direction of the inducing field in Cartesian coordinates
    field_dir = mat_utils.dip_azimuth2cartesian(field_inc, field_dec)

    # Multiply susceptibility model to obtain the x, y, z components of the
    # effective susceptibility contribution from induced magnetization.
    susc_model = np.outer(crust_susc, field_dir)

    rem_model = (rem_susc * mat_utils.dip_azimuth2cartesian(rem_inc, rem_dec))
    
    # Define effective susceptibility model as a vector np.r_[chi_x, chi_y, chi_z]
    plotting_model = susc_model + rem_model
    plotting_model = np.sqrt(np.sum(plotting_model, axis=1) ** 2)

    # populate zmesh
    conversion_factor = 79617.8 # convert to ucgs
    ZZ = np.ones(XX.shape) * plotting_model * conversion_factor

    # insert stripe
    stripe_max = ((ymax+ymin)/2) + stripe_width/2
    stripe_min = ((ymax+ymin)/2) - stripe_width/2
    ZZ[YY>stripe_max] = -plotting_model * conversion_factor
    ZZ[YY<stripe_min] = -plotting_model * conversion_factor

    # output the file
    outputfilename = 'basement_susceptibility_stripe_%.1f.csv' % rem_susc
    hdr = 'x-lims %d %d, y-lims %d %d, cellsize %d, stripe_width %d,crust-susceptibility %.2f, field-inc %.2f, field-dec %.2f, rem-susceptibility %.2f, rem-inc %.2f, rem-dec %.2f' % (xmin,xmax,ymin,ymax,cellsize,stripe_width,crust_susc,field_inc,field_dec,rem_susc,rem_inc,rem_dec)
    np.savetxt(outputfilename,ZZ,delimiter=",",header=hdr,fmt='%.3f')

def build_effective_susc_stripe3(xmin,xmax,ymin,ymax,stripe_width,cellsize,crust_susc,rem_inc,rem_dec,rem_susc,field_inc,field_dec):

    import numpy as np
    from discretize.utils import mkvc, refine_tree_xyz
    from SimPEG.utils import plot2Ddata, model_builder, surface2ind_topo, mat_utils

     # create x-y mesh
    xvec = np.linspace(xmin,xmax,int((xmax+cellsize-xmin)/cellsize))
    yvec = np.linspace(ymin,ymax,int((ymax+cellsize-ymin)/cellsize))
    XX,YY = np.meshgrid(xvec,yvec)

    # Compute the unit direction of the inducing field in Cartesian coordinates
    field_dir = mat_utils.dip_azimuth2cartesian(field_inc, field_dec)

    # Multiply susceptibility model to obtain the x, y, z components of the
    # effective susceptibility contribution from induced magnetization.
    susc_model = np.outer(crust_susc, field_dir)

    rem_model = (rem_susc * mat_utils.dip_azimuth2cartesian(rem_inc, rem_dec))
    
    # Define effective susceptibility model as a vector np.r_[chi_x, chi_y, chi_z]
    plotting_model = susc_model + rem_model
    plotting_model = np.sqrt(np.sum(plotting_model, axis=1) ** 2)

    # populate zmesh
    conversion_factor = 79617.8 # convert to ucgs
    ZZ = np.ones(XX.shape) * plotting_model * conversion_factor

    # insert stripe
    stripe_max = ((ymax+ymin)/2) + stripe_width/2
    stripe_min = ((ymax+ymin)/2) - stripe_width/2
    ZZ[YY>stripe_max] = -plotting_model * conversion_factor
    ZZ[YY<stripe_min] = -plotting_model * conversion_factor

    # insert outer stripes
    stripe1_max = stripe_max + stripe_width
    stripe1_min = stripe_min - stripe_width
    ZZ[YY>stripe1_max] = plotting_model * conversion_factor
    ZZ[YY<stripe1_min] = plotting_model * conversion_factor

    # output the file
    outputfilename = 'basement_susceptibility_stripe3_%.1f.csv' % rem_susc
    hdr = 'x-lims %d %d, y-lims %d %d, cellsize %d, stripe_width %d,crust-susceptibility %.2f, field-inc %.2f, field-dec %.2f, rem-susceptibility %.2f, rem-inc %.2f, rem-dec %.2f' % (xmin,xmax,ymin,ymax,cellsize,stripe_width,crust_susc,field_inc,field_dec,rem_susc,rem_inc,rem_dec)
    np.savetxt(outputfilename,ZZ,delimiter=",",header=hdr,fmt='%.3f')