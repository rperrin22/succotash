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
    outputfilename = 'basement_surface_%d_degrees.csv' % horizon_dip
    hdr = 'x-lims %d %d, y-lims %d %d, cellsize %d, xhinge %d, dip %d, z_base %d' % (xmin,xmax,ymin,ymax,cellsize,xhinge,horizon_dip,z_base)
    np.savetxt(outputfilename,ZZ,delimiter=",",header=hdr)